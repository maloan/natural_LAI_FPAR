# Interpretation:
# Δ(yearmax − yearmean) quantifies whether long-term trends in the
# masked-out (removed) regions are dominated by changes in seasonal
# peak canopy state rather than changes in the annual mean canopy state.
# Positive values indicate relatively stronger peak-season trends.

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
  f_mask_data <- file.path(eval_dir, sprintf("%s_%s_0p25.nc", VAR, "yearmean"))

  # trend slopes for yearmean and yearmax (georef, unmasked)
  f_geo_slope_yearmean <- file.path(geo_dir,
                                    sprintf("%s_georef_yearmean_trend_slope_peryear_0p25.nc", VAR))
  f_geo_slope_yearmax <- file.path(geo_dir,
                                   sprintf("%s_georef_yearmax_trend_slope_peryear_0p25.nc", VAR))
  files <- c(f_geo_data, f_mask_data, f_geo_slope_yearmean, f_geo_slope_yearmax)
  if (!all(file.exists(files))) {
    warning("Missing inputs — skipping ", VAR, " ", MASK, " ", TAU)
    return(invisible(NULL))
  }

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

  # ----------------------------------------------------------------------
  # Baseline mean (1982–2000) for relative normalisation
  # ----------------------------------------------------------------------
  years <- 1982:(1982 + nlyr(r_geo) - 1L)
  idx   <- years >= 1982 & years <= 2000
  if (!any(idx)) {
    warning("Baseline period not covered — skipping ", VAR, " ", MASK, " ", TAU)
    return(invisible(NULL))
  }

  r_base <- mean(r_geo[[idx]], na.rm = TRUE)

  # ----------------------------------------------------------------------
  # Convert slopes to relative (% yr⁻¹)
  # ----------------------------------------------------------------------
  EPS <- 1e-8
  r_slope_mean_rel <- ifel(abs(r_base) < EPS, NA_real_,
                           100 * r_slope_mean / r_base)
  r_slope_max_rel  <- ifel(abs(r_base) < EPS, NA_real_,
                           100 * r_slope_max  / r_base)

  # ----------------------------------------------------------------------
  # Δ-trend = relative(yearmax) − relative(yearmean)
  # ----------------------------------------------------------------------
  r_delta <- r_slope_max_rel - r_slope_mean_rel

  # ----------------------------------------------------------------------
  # Build mask indicating DROPPED pixels
  # (present in unmasked but NA in masked)
  # ----------------------------------------------------------------------
  dropped_mask <- ifel(!is.na(r_geo[[1]]) &
                         is.na(r_mask[[1]]), 1, NA)

  # Δ-trend in dropped region only
  r_dropped_delta <- mask(r_delta, dropped_mask)

  df_trend <- as.data.frame(r_dropped_delta, xy = TRUE, na.rm = TRUE)
  if (nrow(df_trend) < 10) { warning("Too few dropped pixels")
    return(invisible(NULL)) }

  colnames(df_trend) <- c("lon", "lat", "delta")

  # ----------------------------------------------------------------------
  # 1. Δ-trend map
  # ----------------------------------------------------------------------
  p_map <- ggplot(df_trend, aes(lon, lat)) +
    geom_raster(aes(fill = delta)) +
    geom_sf(
      data = coast,
      color = "black",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    scale_fill_scico(
      palette = "bam",
      name = expression(Delta*" trend (% yr"^{-1}*")"),
      oob = scales::squish
    ) +
    scale_x_continuous(breaks = seq(-180, 180, 60), labels = lab_deg) +
    scale_y_continuous(breaks = seq(-90, 90, 30), labels = lab_deg) +
    coord_equal(expand = FALSE) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = sprintf(
        "%s: Δ-trend (yearmax – yearmean) in masked-out region\n(%s, τ = %s)",
        VAR, MASK, TAU
      )
    ) +
    theme_pub()


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
  global_mean_delta <- with(
    df_trend,
    weighted.mean(delta, w = cos(lat * pi / 180), na.rm = TRUE)
  )

  writeLines(
    sprintf(
      "Global mean Δ-trend (yearmax – yearmean) in dropped region: %.6f %% yr-1",
      global_mean_delta
    ),
    con = file.path(OUTDIR, sprintf("summary_delta_%s.txt", MASK))
  )

  # ----------------------------------------------------------------------
  # 3. Zonal mean Δ-trend
  # ----------------------------------------------------------------------
  df_zonal <- df_trend |>
    mutate(lat_band = round(lat)) |>
    group_by(lat_band) |>
    summarise(
      mean_delta = weighted.mean(delta, cos(lat * pi / 180), na.rm = TRUE),
      .groups = "drop"
    )

  p_zonal <- ggplot(df_zonal, aes(mean_delta, lat_band)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_path(linewidth = 0.7, color = "black") +
    scale_y_continuous(labels = lab_deg) +
    labs(
      x = "Δ-trend (yearmax – yearmean) (% yr⁻¹)",
      y = "Latitude (°)",
      title = sprintf("%s: Zonal Δ-trend in masked-out region (%s)", VAR, MASK)
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
