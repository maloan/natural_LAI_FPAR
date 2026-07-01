# =============================================================================
# weighted_means.R — Helper functions for computing weighted means and
# statistics on raster data.
# =============================================================================

global_wmean_series <- function(r, area, year0 = 1982L) {
  # Compute a time series of area-weighted means from a multi-layer raster,
  # using an area raster for weighting.
  terra::compareGeom(r, area, stopOnError = TRUE)
  yrs <- year0:(year0 + terra::nlyr(r) - 1L)
  vals <- vector("numeric", length = terra::nlyr(r))

  for (i in seq_len(terra::nlyr(r))) {
    x <- r[[i]]
    ok <- is.finite(x) & is.finite(area) & area > 0
    num <- as.numeric(terra::global(terra::ifel(ok, x * area, NA), "sum", na.rm = TRUE)[1, 1])
    den <- as.numeric(terra::global(terra::ifel(ok, area, NA), "sum", na.rm = TRUE)[1, 1])
    vals[i] <- ifelse(is.finite(den) &&
      den > 0, num / den, NA_real_)
  }

  data.frame(year = yrs, value = vals)
}

zonal_wmean_latbands <- function(r,
                                 area,
                                 band_deg = 1L,
                                 scale_factor = 1) {
  # Compute zonal weighted means for latitude bands from a raster, using an area
  # raster for weighting. The latitude bands are defined by the specified band
  # width in degrees.
  if (terra::nlyr(r) > 1) {
    r <- r[[1]]
  }

  terra::compareGeom(r, area, stopOnError = TRUE)
  lat <- terra::init(area, "y")
  zone <- floor((lat + 90) / band_deg) + 1

  ok <- is.finite(r) & is.finite(area) & area > 0
  num <- terra::ifel(ok, r * area, NA)
  den <- terra::ifel(ok, area, NA)

  z_num <- terra::zonal(num, zone, "sum", na.rm = TRUE)
  names(z_num) <- c("zone", "num_sum")
  z_den <- terra::zonal(den, zone, "sum", na.rm = TRUE)
  names(z_den) <- c("zone", "den_sum")

  z <- merge(z_num, z_den, by = "zone", all = TRUE)
  data.frame(
    lat_band = -90 + (z$zone - 0.5) * band_deg,
    value = ifelse(
      is.finite(z$num_sum) & is.finite(z$den_sum) & z$den_sum > 0,
      scale_factor * (z$num_sum / z$den_sum),
      NA_real_
    ),
    area_km2 = ifelse(is.finite(z$den_sum), z$den_sum, NA_real_)
  )
}

wmean_series <- function(r, w) {
  # Compute a time series of weighted means from a multi-layer raster, using a
  # weight raster for weighting. The weight raster must have the same geometry
  # as the input raster.
  out <- numeric(nlyr(r))
  for (i in seq_len(nlyr(r))) {
    x <- r[[i]]
    ok <- is.finite(x) & is.finite(w) & w > 0
    num <- as.numeric(global(ifel(ok, x * w, NA), "sum", na.rm = TRUE)[1, 1])
    den <- as.numeric(global(ifel(ok, w, NA), "sum", na.rm = TRUE)[1, 1])
    out[i] <- ifelse(is.finite(den) &&
      den > 0, num / den, NA_real_)
  }
  out
}
weighted_stats <- function(x, w) {
  # Compute weighted mean, standard deviation, and total weight (area) from
  # vectors of values and weights. Returns a list with mean, sd, and area.
  dx <- as.data.frame(c(x, w), na.rm = FALSE)
  names(dx) <- c("x", "w")
  dx <- dx |> filter(is.finite(x), is.finite(w), w > 0)
  if (nrow(dx) == 0) {
    return(list(
      mean = NA_real_,
      sd = NA_real_,
      area = NA_real_
    ))
  }
  mu <- weighted.mean(dx$x, dx$w)
  vv <- sum(dx$w * (dx$x - mu)^2) / sum(dx$w)
  list(
    mean = mu,
    sd = sqrt(vv),
    area = sum(dx$w)
  )
}
