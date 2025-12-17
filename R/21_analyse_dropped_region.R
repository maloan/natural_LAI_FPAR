analyse_dropped_region <- function(VAR,
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

  f_geo_data  <- file.path(geo_dir, sprintf("%s_georef_yearmean_0p25.nc", VAR))
  f_mask_data <- file.path(eval_dir, sprintf("%s_yearmean_0p25.nc", VAR))

  # trend slopes for yearmean and yearmax (georef, unmasked)
  f_geo_slope_yearmean <- file.path(geo_dir,
                                    sprintf("%s_georef_yearmean_trend_slope_peryear_0p25.nc", VAR))
  f_geo_slope_yearmax <- file.path(geo_dir,
                                   sprintf("%s_georef_yearmax_trend_slope_peryear_0p25.nc", VAR))

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

  r_slope_mean <- rast(f_geo_slope_yearmean)
  r_slope_max  <- rast(f_geo_slope_yearmax)

  # Δ-trend (global, unmasked): yearmax – yearmean
  r_delta <- r_slope_max - r_slope_mean

  years <- 1982:(1982 + nlyr(r_geo) - 1)

  # ----------------------------------------------------------------------
  # Build mask indicating DROPPED pixels
  # (present in unmasked but NA in masked)
  # ----------------------------------------------------------------------
  dropped_mask <- ifel(!is.na(r_geo[[1]]) &
                         is.na(r_mask[[1]]), 1, NA)

  # Δ-trend in dropped region only
  r_dropped_delta <- mask(r_delta, dropped_mask)

  df_trend <- as.data.frame(r_dropped_delta, xy = TRUE, na.rm = TRUE)
  colnames(df_trend) <- c("lon", "lat", "delta")

  # ----------------------------------------------------------------------
  # 1. Δ-trend map
  # ----------------------------------------------------------------------
  sdev  <- sd(df_trend$delta, na.rm = TRUE)
  # clamp <- SD_K * sdev
  # df_trend$delta_clamped <- pmax(pmin(df_trend$delta, clamp), -clamp)

  p_map <- ggplot(df_trend, aes(lon, lat)) +
    # geom_tile(aes(fill = delta_clamped)) +
    geom_sf(
      data = coast,
      color = "black",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    scale_fill_scico(
      palette = "bam",
      name = expression(Delta == yearmax - yearmean),
      # limits = c(-clamp, clamp),
      oob = scales::squish
    ) +
    scale_x_continuous(breaks = seq(-180, 180, 60), labels = lab_deg) +
    scale_y_continuous(breaks = seq(-90, 90, 30), labels = lab_deg) +
    coord_sf(expand = FALSE) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = sprintf(
        "%s: Δ-trend (yearmax – yearmean) in MASKED-OUT region \n(%s, τ = %s) ",
        VAR,
        MASK,
        TAU
      ),
    ) + theme_pub()

  ggsave(
    file.path(OUTDIR, sprintf(
      "delta_trend_map_dropped_%s.png", MASK
    )),
    p_map,
    width = 7.2,
    height = 3.8,
    dpi = 330
  )

  # ----------------------------------------------------------------------
  # 2. Global mean Δ-trend
  # ----------------------------------------------------------------------
  global_mean_delta <- mean(df_trend$delta, na.rm = TRUE)
  writeLines(
    sprintf(
      "Global mean Δ-trend (yearmax – yearmean) in dropped region: %.6f per year",
      global_mean_delta
    ),
    con = file.path(OUTDIR, sprintf("summary_delta_%s.txt", MASK))
  )

  # ----------------------------------------------------------------------
  # 3. Zonal mean Δ-trend
  # ----------------------------------------------------------------------
  df_zonal <- df_trend |>
    mutate(lat_band = floor(lat)) |>
    group_by(lat_band) |>
    summarise(mean_delta = mean(delta, na.rm = TRUE),
              .groups = "drop")

  p_zonal <- ggplot(df_zonal, aes(mean_delta, lat_band)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_path(linewidth = 0.7, color = "black") +
    scale_y_continuous(labels = lab_deg) +
    labs(
      x = "Δ-trend (yearmax – yearmean) (per year)",
      y = "Latitude (°)",
      title = sprintf("%s: Zonal Δ-trend in masked-out region (%s)", VAR, MASK),
    ) + theme_pub()

  ggsave(
    file.path(OUTDIR, sprintf(
      "delta_trend_zonal_dropped_%s.png", MASK
    )),
    p_zonal,
    width = 5.3,
    height = 4.8,
    dpi = 330
  )

  message("\nFinished dropped-region Δ-trend analysis: ", OUTDIR)
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------
for (tau in c("tau_0.2", "tau_0.1", "tau_0.05")) {
  for (var in c("LAI", "FPAR")) {
    for (mask in c("CCI", "GLC")) {
      analyse_dropped_region(var, mask, TAU = tau)
    }
  }
}
