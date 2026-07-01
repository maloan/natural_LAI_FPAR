# ==============================================================================
# netcdf.R - Helper functions for reading and processing NetCDF files, including
# variable selection, time extraction, and raster conversion.
# ==============================================================================

`%||%` <- function(a, b)
  if (is.null(a))
    b
else
  a

pick_varname <- function(nc, cfg_vars, var) {
  # Determine the appropriate variable name in a NetCDF file based on
  # configuration and common naming conventions.
  vars <- names(nc$var)
  primary <- cfg_vars$nc_var_name_primary %||% var
  fallback <- cfg_vars$nc_var_name_fallback %||% "auto_first_variable"

  var_u <- toupper(var)
  extras <- switch(
    var_u,
    LAI  = c("lai", "LAI", "Leaf_Area_Index"),
    FPAR = c("FPAR", "fPAR", "fpar", "Fpar"),
    character(0)
  )

  candidates <- unique(c(primary, fallback, extras))
  for (nm in candidates) {
    if (identical(nm, "auto_first_variable"))
      next
    if (nm %in% vars) {
      return(nm)
    }
  }

  coord_names <- tolower(
    c(
      "lon",
      "longitude",
      "x",
      "lat",
      "latitude",
      "y",
      "time",
      "time_counter",
      "t",
      "month"
    )
  )
  ok <- Filter(function(nm) {
    v <- nc$var[[nm]]
    if (is.null(v)) {
      return(FALSE)
    }
    if (tolower(nm) %in% coord_names) {
      return(FALSE)
    }
    nd <- length(v$dim)
    nd >= 2
  }, vars)

  if (length(ok)) {
    return(ok[[1]])
  }

  vars[[1L]]
}

get_lonlat <- function(nc, vcfg) {
  # Extract longitude and latitude values from a NetCDF file, using variable
  # names specified in the configuration or common defaults.
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for get_lonlat().", call. = FALSE)
  }

  loncand <- unique(c(vcfg$nc_lon_name, "lon", "longitude", "LONGITUDE"))
  latcand <- unique(c(vcfg$nc_lat_name, "lat", "latitude", "LATITUDE"))

  lonv <- latv <- NULL

  for (d in nc$dim) {
    dn <- tolower(d$name)
    if (is.null(lonv) && dn %in% tolower(loncand)) {
      lonv <- ncdf4::ncvar_get(nc, d$name)
    }
    if (is.null(latv) && dn %in% tolower(latcand)) {
      latv <- ncdf4::ncvar_get(nc, d$name)
    }
  }

  if (is.null(lonv)) {
    for (nm in loncand) {
      if (!is.null(nc$var[[nm]])) {
        lonv <- ncdf4::ncvar_get(nc, nm)
        break
      }
    }
  }
  if (is.null(latv)) {
    for (nm in latcand) {
      if (!is.null(nc$var[[nm]])) {
        latv <- ncdf4::ncvar_get(nc, nm)
        break
      }
    }
  }

  list(lon = lonv, lat = latv)
}



same_grid <- function(x, y) {
  # Check if two SpatRaster objects have the same grid (CRS, resolution, extent)
  tryCatch({
    identical(terra::crs(x), terra::crs(y)) &&
      isTRUE(all.equal(terra::res(x), terra::res(y))) &&
      isTRUE(all.equal(terra::ext(x), terra::ext(y)))
  }, error = function(e)
    FALSE)
}

transpose_lonlat <- function(arr, lon_len, lat_len) {
  # Transpose a 2D array if its dimensions do not match the expected lon/lat
  # lengths
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
  # Rotate a SpatRaster if its longitude values exceed 180 degrees (i.e., 0-360
  # range)
  if (!is.null(lon_vals) && length(lon_vals) > 0L) {
    if (suppressWarnings(max(lon_vals, na.rm = TRUE)) > 180) {
      return(terra::rotate(r))
    }
    return(r)
  }

  ex <- terra::ext(r)
  if (is.finite(ex[2]) &&
      ex[2] > 180 && is.finite(ex[1]) && ex[1] >= 0) {
    return(terra::rotate(r))
  }
  r
}

.find_time_axis <- function(nc, varnm) {
  # Identify the time axis in a NetCDF file for a given variable, returning its
  # position and corresponding dates.
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for NetCDF time handling.",
         call. = FALSE)
  }

  dims <- nc$var[[varnm]]$dim
  dnames <- tolower(vapply(dims, `[[`, character(1), "name"))

  tpos <- which(dnames %in% c("time", "time_counter", "t", "month"))

  if (!length(tpos)) {
    for (i in seq_along(dims)) {
      u <- try(ncdf4::ncatt_get(nc, dims[[i]]$name, "units")$value,
               silent = TRUE)
      if (is.character(u) &&
          grepl("since", u, ignore.case = TRUE)) {
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
  # Extract time values from a NetCDF file following CF conventions, returning
  # POSIXct dates.
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for NetCDF time handling.",
         call. = FALSE)
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

  m <- regexec("^\\s*(seconds|minutes|hours|days)\\s+since\\s+(.+)$",
               tolower(units))
  mm <- regmatches(tolower(units), m)[[1]]
  if (length(mm) < 3) {
    return(NULL)
  }

  step <- mm[2]
  origin_str <- mm[3]

  origin <- as.POSIXct(if (grepl(":", origin_str)) {
    origin_str
  } else {
    paste0(origin_str, " 00:00:00")
  }, tz = "UTC")

  mult <- switch(
    step,
    seconds = 1,
    minutes = 60,
    hours = 3600,
    days = 86400,
    86400
  )

  as.POSIXct(as.numeric(origin) + tv * mult,
             origin = "1970-01-01",
             tz = "UTC")
}

extract_time_slice <- function(arr, varnm, nc, ym_fallback) {
  # Extract a time slice from a NetCDF variable array based on a year-month
  # fallback string (YYYYMM). If no time axis is found, return the full array.
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
    message(
      sprintf(
        "[extract_time_slice]\n      Using first time slice for %s (no CF match for %s).",
        varnm,
        ym_fallback
      )
    )
    return(do.call(`[`, c(list(arr), ind, list(drop = TRUE))))
  }

  src <- try(nc$filename, silent = TRUE)
  message(
    sprintf(
      "[extract_time_slice]\n    No time dim found; using full array (ym=%s, src=%s, var=%s).",
      ym_fallback,
      if (inherits(src, "try-error"))
        "<unknown>"
      else
        src,
      varnm
    )
  )
  arr
}


agg005_to_025_aw <- function(r005, area005, ref025) {
  # Aggregate a 0.05-degree raster to a 0.25-degree raster using area-weighted
  # averaging, ensuring the grids are nested and aligned.
  stopifnot(terra::nlyr(r005) == 1)

  if (!same_grid(r005, area005)) {
    area005 <- terra::cellSize(r005, unit = "km")
  }

  r5 <- terra::res(r005)
  r25 <- terra::res(ref025)

  fx <- r25[1] / r5[1]
  fy <- r25[2] / r5[2]
  if (!all(abs(c(fx, fy) - round(c(fx, fy))) < 1e-9)) {
    stop(sprintf("Grids are not nested: ref025/res(r005) = (%.8f, %.8f)", fx, fy))
  }

  fact <- c(as.integer(round(fx)), as.integer(round(fy)))

  w <- terra::ifel(is.finite(r005), area005, 0)
  num <- terra::aggregate(r005 * w,
                          fact = fact,
                          fun = sum,
                          na.rm = TRUE)
  den <- terra::aggregate(w,
                          fact = fact,
                          fun = sum,
                          na.rm = TRUE)
  out <- terra::ifel(den == 0, NA, num / den)

  if (!same_grid(out, ref025)) {
    out <- terra::resample(out, ref025, method = "bilinear")
  }
  out
}

nc_month_to_raster <- function(nc_file,
                               ym,
                               var,
                               vcfg,
                               ref,
                               method = "bilinear",
                               extent_global = terra::ext(-180, 180, -90, 90),
                               crs_out = terra::crs(ref),
                               strict_time = TRUE) {
  # Read a specific month from a NetCDF file and return it as a SpatRaster,
  # aligning it to a reference raster.
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for nc_month_to_raster().",
         call. = FALSE)
  }

  stopifnot(is.character(nc_file),
            length(nc_file) == 1L,
            file.exists(nc_file))
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
    stop("Could not determine lon/lat coordinates in NetCDF file: ",
         nc_file,
         call. = FALSE)
  }

  if (length(dim(arr)) == 2L) {
    arr <- transpose_lonlat(arr, length(ll$lon), length(ll$lat))
    r <- terra::rast(
      arr,
      xmin = min(ll$lon),
      xmax = max(ll$lon),
      ymin = min(ll$lat),
      ymax = max(ll$lat),
      crs = crs_out
    )
  } else {
    stop("Expected a 2D slice after time extraction for file: ",
         nc_file,
         call. = FALSE)
  }

  r <- rotate_if_360(r, lon_vals = ll$lon)
  if (!terra::compareGeom(r, ref, stopOnError = FALSE)) {
    r <- align_to_template(r, ref, method = method)
  }
  r
}
