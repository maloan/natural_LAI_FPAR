## =============================================================================
## geom.R — Core geometry + NetCDF monthly slice → aligned SpatRaster
## =============================================================================

source(here("R", "helpers", "utils.R"))

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
  # Prefer actual coordinate values when available
  if (!is.null(lon_vals) && length(lon_vals) > 0L) {
    if (suppressWarnings(max(lon_vals, na.rm = TRUE)) > 180) {
      return(terra::rotate(r))
    }
    return(r)
  }

  # Fallback heuristic: raster extent in [0, 360]
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
  # returns: list(tpos = integer index in arr dims, dates = POSIXct or NULL)
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for NetCDF time handling.", call. = FALSE)
  }

  dims <- nc$var[[varnm]]$dim
  dnames <- tolower(vapply(dims, `[[`, character(1), "name"))

  # 1) direct name match
  tpos <- which(dnames %in% c("time", "time_counter", "t", "month"))

  # 2) CF-units "since" on dim variable
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

  # ensure time component
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

  # If we have CF dates, pick closest to requested month
  if (length(tpos) && !is.null(dates)) {
    ym_target <- suppressWarnings(as.Date(paste0(ym_fallback, "01"), "%Y%m%d"))
    if (is.finite(ym_target)) {
      idx <- which.min(abs(as.Date(dates) - ym_target))
      ind <- rep(list(TRUE), length(d))
      ind[[tpos]] <- idx
      return(do.call(`[`, c(list(arr), ind, list(drop = TRUE))))
    }
  }

  # If we have a time axis but no dates, take first slice
  if (length(tpos)) {
    ind <- rep(list(TRUE), length(d))
    ind[[tpos]] <- 1L
    message(sprintf(
      "[extract_time_slice]
      Using first time slice for %s (no CF match for %s).", varnm, ym_fallback
    ))
    return(do.call(`[`, c(list(arr), ind, list(drop = TRUE))))
  }

  # No time axis detected: return arr unchanged (forgiving)
  src <- try(nc$filename, silent = TRUE)
  message(sprintf(
    "[extract_time_slice]
    No time dim found; using full array (ym=%s, src=%s, var=%s).",
    ym_fallback,
    if (inherits(src, "try-error")) "<unknown>" else src,
    varnm
  ))
  arr
}

# ------------------------------------------------------------------------------
# Area-weighted aggregation 0.05° → 0.25° (nested-grid assumption)
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
# NetCDF monthly slice → aligned SpatRaster
#   - pick_varname(), get_lonlat(), align_to() are expected from utils.R
# ------------------------------------------------------------------------------

nc_month_to_raster <- function(nc_file,
                               ym,
                               VAR,
                               Vcfg,
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

  VAR <- toupper(VAR)

  nc <- ncdf4::nc_open(nc_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  varnm <- pick_varname(nc, Vcfg, VAR)

  arr <- ncdf4::ncvar_get(nc, varnm)
  arr <- extract_time_slice(arr, varnm, nc, ym)

  if (is.null(arr)) {
    msg <- sprintf(
      "Time slice '%s' missing in '%s' (var=%s).", ym, basename(nc_file), varnm
    )
    if (isTRUE(strict_time)) stop(msg, call. = FALSE) else return(NULL)
  }

  ll <- get_lonlat(nc, Vcfg)

  if (length(dim(arr)) == 2L && !is.null(ll$lon) && !is.null(ll$lat)) {
    arr <- transpose_lonlat(arr, length(ll$lon), length(ll$lat))
  }

  atts <- ncdf4::ncatt_get(nc, varnm)
  fillv <- atts$`_FillValue` %||% atts$missing_value %||% NA
  if (is.finite(fillv)) arr[arr == fillv] <- NA_real_

  r <- terra::rast(arr)
  terra::ext(r) <- extent_global
  terra::crs(r) <- crs_out

  r <- rotate_if_360(r, lon_vals = ll$lon)
  align_to(r, ref, method = method)
}
