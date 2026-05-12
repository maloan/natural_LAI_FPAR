## =============================================================================
## weighted_means.R — Shared weighted mean helpers
## =============================================================================

global_wmean_series <- function(r, area, year0 = 1982L) {
  terra::compareGeom(r, area, stopOnError = TRUE)
  yrs <- year0:(year0 + terra::nlyr(r) - 1L)
  vals <- vector("numeric", length = terra::nlyr(r))

  for (i in seq_len(terra::nlyr(r))) {
    x <- r[[i]]
    ok <- is.finite(x) & is.finite(area) & area > 0
    num <- as.numeric(terra::global(terra::ifel(ok, x * area, NA), "sum", na.rm = TRUE)[1, 1])
    den <- as.numeric(terra::global(terra::ifel(ok, area, NA), "sum", na.rm = TRUE)[1, 1])
    vals[i] <- ifelse(is.finite(den) && den > 0, num / den, NA_real_)
  }

  data.frame(year = yrs, value = vals)
}

zonal_wmean_latbands <- function(r, area, band_deg = 1L, scale_factor = 1) {
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
