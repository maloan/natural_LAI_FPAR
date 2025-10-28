## =============================================================================
# geom.R — Core geometry and NetCDF handling utilities
#
# Purpose
#   Shared low-level tools for raster alignment, CF-compliant time extraction,
#   and area-weighted aggregation used across preprocessing and evaluation steps.
#
# Functions
#   %||%()                    — Null-coalescing operator (x if not NULL else y)
#   same_grid()               — Test CRS, extent, and resolution equality
#   maybe_transpose_lonlat()  — Swap lon/lat array axes if order mismatched
#   maybe_flip_lat()          — Flip array vertically if latitudes descend
#   rotate_if_360()           — Convert 0–360° → −180–180° longitude convention
#   .extract_cf_times()       — Parse CF-compliant time variable from NetCDF
#   extract_time_slice()      — Extract nearest CF-time slice by YYYYMM fallback
#   agg005_to_025_aw()        — Area-weighted 0.05°→0.25° aggregation
#   stack_aligned()           — Load and align rasters to a reference grid
#   get_area_for()            — Cached per-grid area computation
#
# Inputs
#   Raster or NetCDF objects from individual workflow steps.
#
# Outputs
#   terra objects, numeric arrays, or POSIXct time vectors.
#
# Dependencies
#   terra, ncdf4
#
# Processing overview
#   1) Provide robust raster alignment and dimension handling.
#   2) Implement efficient CF-time parsing and fallback logic.
#   3) Support area-weighted upscaling and memoized area caching.
## =============================================================================

same_grid <- function(x, y) {
  tryCatch({
    (terra::crs(x) == terra::crs(y)) &&
      isTRUE(all.equal(terra::res(x), terra::res(y))) &&
      isTRUE(all.equal(terra::ext(x), terra::ext(y)))
  }, error = function(e)
    FALSE)
}

maybe_transpose_lonlat <- function(arr, lon_len, lat_len) {
  d <- dim(arr)
  if (length(d) != 2L)
    return(arr)
  if (d[1] == lon_len && d[2] == lat_len)
    return(arr)
  if (d[1] == lat_len &&
      d[2] == lon_len)
    return(aperm(arr, c(2, 1)))
  arr
}

maybe_flip_lat <- function(arr, lat_vals) {
  if (is.null(lat_vals) || length(dim(arr)) != 2L)
    return(arr)
  if (length(lat_vals) > 1 && all(diff(lat_vals) < 0))
  {
    arr[, ncol(arr):1, drop = FALSE]
  }
  else {
    arr
  }
}

rotate_if_360 <- function(r, lon_vals = NULL) {
  try({
    if (!is.null(lon_vals) && max(lon_vals, na.rm = TRUE) > 180)
    {
      return(terra::rotate(r))
    }
  }, silent = TRUE)
  r
}

.extract_cf_times <- function(nc) {
  all_vars <- names(nc$var)
  cand <- intersect(all_vars, c("time", "time_counter", "t", "month"))
  if (!length(cand)) {
    cand <- Filter(function(v) {
      at <- try(ncdf4::ncatt_get(nc, v, "units")$value, silent = TRUE)
      is.character(at) && grepl("since", at, ignore.case = TRUE)
    }, all_vars)
  }
  if (!length(cand))
    return(NULL)
  tvar <- cand[[1]]
  tv   <- as.numeric(ncdf4::ncvar_get(nc, tvar))
  units <- try(ncdf4::ncatt_get(nc, tvar, "units")$value, silent = TRUE)
  if (inherits(units, "try-error") ||
      !is.character(units))
    return(NULL)

  m <- regexec("^\\s*(seconds|minutes|hours|days)\\s+since\\s+([^\\s]+)",
               tolower(units))
  mm <- regmatches(tolower(units), m)[[1]]
  if (length(mm) < 3)
    return(NULL)

  step   <- mm[2]
  origin <- as.POSIXct(paste0(mm[3], if (!grepl(":", mm[3]))
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
  if (length(dim(arr)) <= 2L)
    return(arr)

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

  dates <- .extract_cf_times(nc)
  if (length(tpos) && !is.null(dates)) {
    ym_target <- try(as.Date(paste0(ym_fallback, "01"), "%Y%m%d"), silent = TRUE)
    if (!inherits(ym_target, "try-error")) {
      idx <- which.min(abs(as.Date(dates) - ym_target))
      slicer <- rep(list(bquote()), length(dim(arr)))
      slicer[[tpos[1]]] <- idx
      return(do.call(`[`, c(list(arr), slicer, list(drop = TRUE))))
    }
  }

  if (!length(tpos)) {
    src <- try(nc$filename, silent = TRUE)
    message(
      sprintf(
        "[extract_time_slice] No time dim found; using filename token %s (src: %s, var: %s)",
        ym_fallback,
        if (inherits(src, "try-error")) {
          "<unknown>"
        } else {
          src
        },
        varnm
      )
    )
    return(arr)
  }

  message(
    sprintf(
      "[extract_time_slice] Using first time slice for %s (no match for %s).",
      varnm,
      ym_fallback
    )
  )
  tpos <- tpos[1]
  slicer <- rep(list(quote(expr =)), length(dim(arr)))
  slicer[[tpos]] <- 1L
  do.call(`[`, c(list(arr), slicer, list(drop = TRUE)))
}

agg005_to_025_aw <- function(r005, area005, ref025) {
  stopifnot(terra::nlyr(r005) == 1)
  if (!same_grid(r005, area005))
  {
    area005 <- terra::cellSize(r005, unit = "km")
  }

  r5_res  <- terra::res(r005)
  r25_res <- terra::res(ref025)
  fx <- r25_res[1] / r5_res[1]
  fy <- r25_res[2] / r5_res[2]
  if (!all(abs(c(fx, fy) - round(c(fx, fy))) < 1e-9))
  {
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

  if (!same_grid(out, ref025))
  {
    out <- terra::resample(out, ref025, method = "bilinear")
  }

  out
}

stack_aligned <- function(files, template, method = "bilinear") {
  rs <- lapply(files, function(f) {
    r <- terra::rast(f)
    if (!same_grid(r, template))
    {
      r <- terra::resample(r, template, method = method)
    }
    r
  })
  terra::rast(rs)
}

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
