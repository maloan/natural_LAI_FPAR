`%||%` <- function(x, y)
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }

same_grid <- function(x, y) {
  tryCatch({
    (terra::crs(x) == terra::crs(y)) &&
      isTRUE(all.equal(terra::res(x), terra::res(y))) &&
      isTRUE(all.equal(terra::ext(x), terra::ext(y)))
  }, error = function(e)
    FALSE)
}

.add_overlays <- function(r, lwd = 1.2, col = "black") {
  if (!exists("cfg") || is.null(cfg$aois)) return()
  for (aoi in cfg$aois) {
    if (inherits(aoi, "SpatVector") || inherits(aoi, "sf")) {
      lines(vect(aoi), lwd = lwd, col = col)
    }
  }
}

maybe_transpose_lonlat <- function(arr, lon_len, lat_len) {
  d <- dim(arr)
  if (length(d) != 2L)
    return(arr)
  # guess whether arr is [lon,lat] or [lat,lon]
  if (d[1] == lon_len &&
      d[2] == lat_len)
    return(arr)       # [lon,lat]
  if (d[1] == lat_len && d[2] == lon_len)
    return(aperm(arr, c(2, 1)))
  arr  # unknown, leave as-is
}

maybe_flip_lat <- function(arr, lat_vals) {
  if (is.null(lat_vals) || length(dim(arr)) != 2L)
    return(arr)
  if (length(lat_vals) > 1 && all(diff(lat_vals) < 0)) {
    # flip top<->bottom so latitude increases from bottom to top
    arr[, ncol(arr):1, drop = FALSE]
  } else {
    arr
  }
}


rotate_if_360 <- function(r, lon_vals = NULL) {
  try({
    if (!is.null(lon_vals) &&
        max(lon_vals, na.rm = TRUE) > 180)
      return(terra::rotate(r))
  }, silent = TRUE)
  r
}

.extract_cf_times <- function(nc) {
  # find a CF-like time variable
  all_vars <- names(nc$var)
  cand <- intersect(all_vars, c("time", "time_counter", "t", "month"))
  if (!length(cand)) {
    # look for any var with CF "since"
    cand <- Filter(function(v) {
      at <- try(ncdf4::ncatt_get(nc, v, "units")$value, silent = TRUE)
      is.character(at) && grepl("since", at, ignore.case = TRUE)
    }, all_vars)
  }
  if (!length(cand))
    return(NULL)
  tvar <- cand[[1]]
  tv  <- as.numeric(ncdf4::ncvar_get(nc, tvar))
  units <- try(ncdf4::ncatt_get(nc, tvar, "units")$value, silent = TRUE)
  if (inherits(units, "try-error") ||
      !is.character(units))
    return(NULL)
  # parse "X since YYYY-MM-DD ..."
  m <- regexec("^\\s*(seconds|minutes|hours|days)\\s+since\\s+([^\\s]+)",
               tolower(units))
  mm <- regmatches(tolower(units), m)[[1]]
  if (length(mm) < 3)
    return(NULL)
  step <- mm[2]
  origin <- mm[3]
  origin <- as.POSIXct(paste0(origin, if (!grepl(":", origin))
    " 00:00:00"
    else
      ""), tz = "UTC")
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
  # If data is already 2D, nothing to slice
  if (length(dim(arr)) <= 2L)
    return(arr)
  
  # Which dim is time?
  dims <- nc$var[[varnm]]$dim
  dnames <- tolower(vapply(dims, `[[`, character(1), "name"))
  tpos <- which(dnames %in% c("time", "time_counter", "t", "month"))
  if (!length(tpos)) {
    # heuristic: look for a dim whose var has CF time units
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
  
  # Try to match YYYYMM in filename to the CF time axis
  dates <- .extract_cf_times(nc)
  if (length(tpos) && !is.null(dates)) {
    ym_target <- try(as.Date(paste0(ym_fallback, "01"), "%Y%m%d"), silent = TRUE)
    if (!inherits(ym_target, "try-error")) {
      # choose index with the smallest |date - ym_target|
      idx <- which.min(abs(as.Date(dates) - ym_target))
      slicer <- rep(list(bquote()), length(dim(arr)))
      slicer[[tpos[1]]] <- idx
      return(do.call(`[`, c(list(arr), slicer, list(drop = TRUE))))
    }
  }
  
  # Fallbacks
  if (!length(tpos)) {
    src <- try(nc$filename, silent = TRUE)
    message(
      sprintf(
        "[extract_time_slice] No time dim found; using filename token %s (src: %s, var: %s)",
        ym_fallback,
        if (inherits(src, "try-error"))
          "<unknown>"
        else
          src,
        varnm
      )
    )
    return(arr)
  }
  # Time exists but no parseable match — take first slice, warn
  message(
    sprintf(
      "[extract_time_slice] Using first time slice for %s (could not match %s).",
      varnm,
      ym_fallback
    )
  )
  tpos <- tpos[1]
  nd <- length(dim(arr))
  slicer <- rep(list(quote(expr = )), nd)
  slicer[[tpos]] <- 1L
  do.call(`[`, c(list(arr), slicer, list(drop = TRUE)))
}



# 0.05° -> 0.25° area-weighted aggregation
agg005_to_025_aw <- function(r005, area005, ref025) {
  stopifnot(terra::nlyr(r005) == 1)
  
  if (!same_grid(r005, area005)) {
    area005 <- terra::cellSize(r005, unit = "km")
  }
  
  r5_res  <- terra::res(r005)
  r25_res <- terra::res(ref025)
  fx <- r25_res[1] / r5_res[1]
  fy <- r25_res[2] / r5_res[2]
  if (!all(abs(c(fx, fy) - round(c(fx, fy))) < 1e-9)) {
    stop(sprintf("Grids are not nested: ref025/res(r005) = (%.8f, %.8f)", fx, fy))
  }
  fact <- c(as.integer(round(fx)), as.integer(round(fy)))
  
  w   <- terra::ifel(is.finite(r005), area005, 0)
  num <- terra::aggregate(r005 * w,
                          fact = fact,
                          fun = sum,
                          na.rm = TRUE)
  den <- terra::aggregate(w,
                          fact = fact,
                          fun = sum,
                          na.rm = TRUE)
  out <- num / den
  out <- terra::ifel(den == 0, NA, out)
  
  if (!same_grid(out, ref025)) {
    out <- terra::resample(out, ref025, method = "bilinear")  # preserves block means
  }
  out
}



# Stack aligned to a template
stack_aligned <- function(files,
                          template,
                          method = c("near", "bilinear")[2]) {
  rs <- lapply(files, function(f) {
    r <- terra::rast(f)
    if (!same_grid(r, template)) {
      r <- terra::resample(r, template, method = method)
    }
    r
  })
  terra::rast(rs)
}

# Memoized area for arbitrary grids
.get_area_cache <- new.env(parent = emptyenv())
.sig <- function(r)
  paste0(
    terra::crs(r),
    "|",
    paste(terra::res(r), collapse = ","),
    "|",
    paste(round(
      c(
        terra::xmin(r),
        terra::xmax(r),
        terra::ymin(r),
        terra::ymax(r)
      ), 6
    ), collapse = ",")
  )
get_area_for <- function(r) {
  key <- .sig(r)
  if (exists(key, envir = .get_area_cache, inherits = FALSE))
    return(get(key, envir = .get_area_cache))
  a <- try(terra::cellSize(r, unit = "km"), silent = TRUE)
  if (inherits(a, "try-error")) {
    a <- terra::rast(r)
    terra::values(a) <- NA_real_
  }
  assign(key, a, envir = .get_area_cache)
  a
}
