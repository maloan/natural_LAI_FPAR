
## =============================================================================
# 17_analyse_trend_mask.R
## =============================================================================
analyse_dropped_region <- function(VAR,
                                   METRIC,
                                   MASK,
                                   TAU = "tau_0.1",
                                   OUTDIR = "analysis/trends_dropped",
                                   SD_K = 2) {
  suppressPackageStartupMessages({
    library(terra)
    library(dplyr)
    library(ggplot2)
    library(scico)
    library(sf)
    library(rnaturalearth)
    library(scales)
    library(patchwork)
  })

  # ----------------------------------------------------------------------
  # File paths
  # ----------------------------------------------------------------------
  geo_dir   <- file.path("analysis/unmasked/0p25")
  eval_dir  <- file.path("output", TAU, "eval", sprintf("trend_%s_%s", VAR, MASK))
  f_geo_data  <- file.path(geo_dir, sprintf("%s_georef_%s_0p25.nc", VAR, METRIC))
  f_geo_slope <- file.path(geo_dir,
                           sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", VAR, METRIC))
  f_mask_data  <- file.path(eval_dir, sprintf("%s_%s_0p25.nc", VAR, METRIC))
  f_mask_slope <- file.path(eval_dir,
                            sprintf("%s_%s_trend_slope_peryear_0p25.nc", VAR, METRIC))

  # ----------------------------------------------------------------------
  # Output directory
  # ----------------------------------------------------------------------
  OUTDIR <- file.path(OUTDIR, sprintf("%s/%s", TAU, VAR))
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

  # ----------------------------------------------------------------------
  # Theme & coastlines
  # ----------------------------------------------------------------------
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
  theme_pub <- function() {
    theme_bw(base_size = 12) +
      theme(
        panel.grid.major = element_line(color = "grey87", linewidth = 0.3),
        axis.text  = element_text(size = 9),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 13, face = "bold")
      )
  }
  lab_deg <- scales::label_number(suffix = "°")

  # ----------------------------------------------------------------------
  # Load rasters
  # ----------------------------------------------------------------------
  r_geo   <- rast(f_geo_data)
  r_mask  <- rast(f_mask_data)
  r_slope <- rast(f_geo_slope)
  years <- 1982:(1982 + nlyr(r_geo) - 1)
  # ----------------------------------------------------------------------
  # Build mask indicating DROPPED pixels
  # (present in unmasked but NA in masked)
  # ----------------------------------------------------------------------
  dropped_mask <- ifel(!is.na(r_geo[[1]]) &
                         is.na(r_mask[[1]]), 1, NA)
  # Trend in dropped region only
  r_dropped_trend <- mask(r_slope, dropped_mask)
  df_trend <- as.data.frame(r_dropped_trend, xy = TRUE, na.rm = TRUE)
  colnames(df_trend) <- c("lon", "lat", "slope")
  # ----------------------------------------------------------------------
  # 1. Trend map
  # ----------------------------------------------------------------------
  sdev  <- sd(df_trend$slope, na.rm = TRUE)
  clamp <- SD_K * sdev
  df_trend$slope_clamped <- pmax(pmin(df_trend$slope, clamp), -clamp)
  p_map <- ggplot(df_trend, aes(lon, lat)) +
    geom_tile(aes(fill = slope_clamped)) +
    geom_sf(
      data = coast,
      color = "black",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    scale_fill_scico(
      palette = "bam",
      name = sprintf("%s trend (dropped region)", VAR),
      limits = c(-clamp, clamp),
      oob = scales::squish
    ) +
    scale_x_continuous(breaks = seq(-180, 180, 60), labels = lab_deg) +
    scale_y_continuous(breaks = seq(-90, 90, 30), labels = lab_deg) +
    coord_sf(expand = FALSE) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = sprintf("%s %s: Trend in MASKED-OUT region (%s)", VAR, METRIC, MASK)
    ) + theme_pub()
  ggsave(
    file.path(
      OUTDIR,
      sprintf("trend_map_dropped_%s_%s.png", METRIC, MASK)
    ),
    p_map,
    width = 7.2,
    height = 3.8,
    dpi = 330
  )

  # ----------------------------------------------------------------------
  # 2. Global mean trend
  # ----------------------------------------------------------------------
  global_mean_trend <- mean(df_trend$slope, na.rm = TRUE)
  writeLines(
    sprintf(
      "Global mean trend in dropped region: %.6f per year",
      global_mean_trend
    ),
    con = file.path(OUTDIR, "summary.txt")
  )
  # ----------------------------------------------------------------------
  # 3. Zonal mean trend
  # ----------------------------------------------------------------------
  df_zonal <- df_trend |>
    mutate(lat_band = floor(lat)) |>
    group_by(lat_band) |>
    summarise(mean_trend = mean(slope, na.rm = TRUE),
              .groups = "drop")
  p_zonal <- ggplot(df_zonal, aes(mean_trend, lat_band)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_path(linewidth = 0.7, color = "black") +
    scale_y_continuous(labels = lab_deg) +
    labs(
      x = sprintf("%s trend (per year)", VAR),
      y = "Latitude (°)",
      title = sprintf("%s %s: Zonal trend of masked-out region", VAR, METRIC)
    ) + theme_pub()
  ggsave(
    file.path(
      OUTDIR,
      sprintf("trend_zonal_dropped_%s_%s.png", METRIC, MASK)
    ),
    p_zonal,
    width = 5.3,
    height = 4.8,
    dpi = 330
  )

  # ----------------------------------------------------------------------
  # 4. Time series for the dropped regions
  # ----------------------------------------------------------------------
  # For each year: mask yearly raster -> compute global mean
  dropped_ts <- tibble(year = years, value = sapply(1:nlyr(r_geo), function(i) {
    r_year <- mask(r_geo[[i]], dropped_mask)
    mean(values(r_year), na.rm = TRUE)
  }))

  p_ts <- ggplot(dropped_ts, aes(year, value)) +
    geom_line(linewidth = 0.55) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "red",
      linewidth = 0.55
    ) +
    labs(
      x = "Year",
      y = sprintf("%s (global mean)", VAR),
      title = sprintf("%s %s: Annual series (masked-out region)", VAR, METRIC),
      subtitle = sprintf("%s mask", MASK)
    ) +
    theme_pub()
  ggsave(
    file.path(
      OUTDIR,
      sprintf("timeseries_dropped_%s_%s.png", METRIC, MASK)
    ),
    p_ts,
    width = 6.4,
    height = 4.2,
    dpi = 330
  )
  message("\nFinished dropped-region analysis: ", OUTDIR)
}
# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------
for (tau in c("tau_0.05", "tau_0.1", "tau_0.2")) {
  for (var in c("LAI", "FPAR")) {
    for (met in c("yearmean", "yearmax")) {
      for (mask in c("CCI", "GLC")) {
        analyse_dropped_region(var, met, mask, TAU = tau)
      }
    }
  }
}
