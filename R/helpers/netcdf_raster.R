## =============================================================================
## netcdf_raster.R — Core geometry + NetCDF monthly slice -> aligned SpatRaster
## =============================================================================

source(here("R", "helpers", "paths.R"))
source(here("R", "helpers", "files.R"))
source(here("R", "helpers", "netcdf.R"))

# ------------------------------------------------------------------------------
# Basic grid checks
# ------------------------------------------------------------------------------

same_grid <- function(x, y) {
  tryCatch(
    {
      identical(terra::crs(x), terra::crs(y)) &&
        isTRUE(all.equal(terra::res(x), terra::res(y))) &&
        isTRUE(all.equal(terra::ext(x), terra::ext(y)))
    },
    error = function(e) FALSE
  )
}

transpose_lonlat <- function(arr, lon_len, lat_len) {
  d <- dim(arr)
  if (length(d) != 2L) {
    return(arr)
  }
  if (d[1] == lon_len && d[2] == lat_len) {
    return(arr)
  }
  if (d[1] == lat_len && d[2] == lon_len) {
    return(aperm(arr, c(2, 1)))
  }
  arr
}

rotate_if_360 <- function(r, lon_vals = NULL) {
  if (!is.null(lon_vals) && length(lon_vals) > 0L) {
    if (suppressWarnings(max(lon_vals, na.rm = TRUE)) > 180) {
      return(terra::rotate(r))
    }
    return(r)
  }

  ex <- terra::ext(r)
  if (is.finite(ex[2]) && ex[2] > 180 && is.finite(ex[1]) && ex[1] >= 0) {
    return(terra::rotate(r))
  }
  r
}

# ------------------------------------------------------------------------------
# CF-time parsing helpers
# ------------------------------------------------------------------------------

.find_time_axis <- function(nc, varnm) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for NetCDF time handling.", call. = FALSE)
  }

  dims <- nc$var[[varnm]]$dim
  dnames <- tolower(vapply(dims, `[[`, character(1), "name"))

  tpos <- which(dnames %in% c("time", "time_counter", "t", "month"))

  if (!length(tpos)) {
    for (i in seq_along(dims)) {
      u <- try(
        ncdf4::ncatt_get(nc, dims[[i]]$name, "units")$value,
        silent = TRUE
      )
      if (is.character(u) && grepl("since", u, ignore.case = TRUE)) {
        tpos <- i
        break
      }
    }
  }

  list(tpos = if (length(tpos)) {
    tpos[1]
  } else {
    integer(0)
  }, dates = .extract_cf_times(nc))
}

.extract_cf_times <- function(nc) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for NetCDF time handling.", call. = FALSE)
  }

  all_vars <- names(nc$var)
  cand <- intersect(all_vars, c("time", "time_counter", "t", "month"))

  if (!length(cand)) {
    cand <- Filter(function(v) {
      u <- try(ncdf4::ncatt_get(nc, v, "units")$value, silent = TRUE)
      is.character(u) && grepl("since", u, ignore.case = TRUE)
    }, all_vars)
  }
  if (!length(cand)) {
    return(NULL)
  }

  tvar <- cand[[1]]
  tv <- suppressWarnings(as.numeric(ncdf4::ncvar_get(nc, tvar)))
  units <- try(ncdf4::ncatt_get(nc, tvar, "units")$value, silent = TRUE)
  if (!is.character(units)) {
    return(NULL)
  }

  m <- regexec(
    "^\\s*(seconds|minutes|hours|days)\\s+since\\s+(.+)$", tolower(units)
  )
  mm <- regmatches(tolower(units), m)[[1]]
  if (length(mm) < 3) {
    return(NULL)
  }

  step <- mm[2]
  origin_str <- mm[3]

  origin <- as.POSIXct(
    if (grepl(":", origin_str)) {
      origin_str
    } else {
      paste0(origin_str, " 00:00:00")
    },
    tz = "UTC"
  )

  mult <- switch(step,
    seconds = 1,
    minutes = 60,
    hours = 3600,
    days = 86400,
    86400
  )

  as.POSIXct(as.numeric(origin) + tv * mult, origin = "1970-01-01", tz = "UTC")
}

extract_time_slice <- function(arr, varnm, nc, ym_fallback) {
  d <- dim(arr)
  if (length(d) <= 2L) {
    return(arr)
  }

  ax <- .find_time_axis(nc, varnm)
  tpos <- ax$tpos
  dates <- ax$dates

  if (length(tpos) && !is.null(dates)) {
    ym_target <- suppressWarnings(as.Date(paste0(ym_fallback, "01"), "%Y%m%d"))
    if (is.finite(ym_target)) {
      idx <- which.min(abs(as.Date(dates) - ym_target))
      ind <- rep(list(TRUE), length(d))
      ind[[tpos]] <- idx
      return(do.call(`[`, c(list(arr), ind, list(drop = TRUE))))
    }
  }

  if (length(tpos)) {
    ind <- rep(list(TRUE), length(d))
    ind[[tpos]] <- 1L
    message(sprintf(
      "[extract_time_slice]\n      Using first time slice for %s (no CF match for %s).", varnm, ym_fallback
    ))
    return(do.call(`[`, c(list(arr), ind, list(drop = TRUE))))
  }

  src <- try(nc$filename, silent = TRUE)
  message(sprintf(
    "[extract_time_slice]\n    No time dim found; using full array (ym=%s, src=%s, var=%s).",
    ym_fallback,
    if (inherits(src, "try-error")) "<unknown>" else src,
    varnm
  ))
  arr
}

# ------------------------------------------------------------------------------
# Area-weighted aggregation 0.05° -> 0.25° (nested-grid assumption)
# ------------------------------------------------------------------------------

agg005_to_025_aw <- function(r005, area005, ref025) {
  stopifnot(terra::nlyr(r005) == 1)

  if (!same_grid(r005, area005)) {
    area005 <- terra::cellSize(r005, unit = "km")
  }

  r5 <- terra::res(r005)
  r25 <- terra::res(ref025)

  fx <- r25[1] / r5[1]
  fy <- r25[2] / r5[2]
  if (!all(abs(c(fx, fy) - round(c(fx, fy))) < 1e-9)) {
    stop(sprintf(
      "Grids are not nested: ref025/res(r005) = (%.8f, %.8f)", fx, fy
    ))
  }

  fact <- c(as.integer(round(fx)), as.integer(round(fy)))

  w <- terra::ifel(is.finite(r005), area005, 0)
  num <- terra::aggregate(r005 * w, fact = fact, fun = sum, na.rm = TRUE)
  den <- terra::aggregate(w, fact = fact, fun = sum, na.rm = TRUE)
  out <- terra::ifel(den == 0, NA, num / den)

  if (!same_grid(out, ref025)) {
    out <- terra::resample(out, ref025, method = "bilinear")
  }
  out
}

# ------------------------------------------------------------------------------
# NetCDF monthly slice -> aligned SpatRaster
# ------------------------------------------------------------------------------

nc_month_to_raster <- function(nc_file,
                               ym,
                               var,
                               vcfg,
                               ref,
                               method = "bilinear",
                               extent_global = terra::ext(-180, 180, -90, 90),
                               crs_out = terra::crs(ref),
                               strict_time = TRUE) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for nc_month_to_raster().", call. = FALSE)
  }

  stopifnot(is.character(nc_file), length(nc_file) == 1L, file.exists(nc_file))
  stopifnot(inherits(ref, "SpatRaster"))
  stopifnot(is.character(ym), nchar(ym) >= 6)

  var <- toupper(var)

  nc <- ncdf4::nc_open(nc_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  varnm <- pick_varname(nc, vcfg, var)

  arr <- ncdf4::ncvar_get(nc, varnm)
  arr <- extract_time_slice(arr, varnm, nc, ym)

  if (is.null(arr)) {
    msg <- sprintf("Could not read data from %s (var=%s).", nc_file, varnm)
    if (isTRUE(strict_time)) {
      stop(msg, call. = FALSE)
    }
    warning(msg, call. = FALSE)
    return(terra::rast(ref))
  }

  ll <- get_lonlat(nc, vcfg)
  if (is.null(ll$lon) || is.null(ll$lat)) {
    stop("Could not determine lon/lat coordinates in NetCDF file: ", nc_file, call. = FALSE)
  }

  if (length(dim(arr)) == 2L) {
    arr <- transpose_lonlat(arr, length(ll$lon), length(ll$lat))
    r <- terra::rast(
      arr,
      xmin = min(ll$lon), xmax = max(ll$lon),
      ymin = min(ll$lat), ymax = max(ll$lat),
      crs = crs_out
    )
  } else {
    stop("Expected a 2D slice after time extraction for file: ", nc_file, call. = FALSE)
  }

  r <- rotate_if_360(r, lon_vals = ll$lon)
  if (!terra::compareGeom(r, ref, stopOnError = FALSE)) {
    r <- align_to(r, ref, method = method)
  }
  r
}
