
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
  files <- c(f_geo_data, f_geo_slope, f_mask_data)
  if (!all(file.exists(files))) {
    warning("Missing inputs — skipping ",
            VAR, " ", METRIC, " ", MASK, " ", TAU)
    return(invisible(NULL))
  }

  # ----------------------------------------------------------------------
  # Output directory
  # ----------------------------------------------------------------------
  OUTDIR <- file.path(OUTDIR, TAU, VAR, sprintf("%s_%s", METRIC, MASK))
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
  dropped_mask <- ifel(!is.na(r_geo[[1]]) & is.na(r_mask[[1]]), 1, NA)

  # baseline mean for dropped region
  r_base <- baseline_mean_raster(r_geo, years, 1982, 2000)
  r_dropped_trend <- r_slope * dropped_mask
  EPS <- 1e-8
  r_dropped_trend <- ifel(abs(r_base) < EPS, NA_real_,
                          100 * r_dropped_trend / r_base)


  df_trend <- as.data.frame(r_dropped_trend, xy = TRUE, na.rm = TRUE)
  if (nrow(df_trend) < 10) {
    warning("Too few dropped pixels — skipping plots")
    return(invisible(NULL))
  }

  colnames(df_trend) <- c("lon", "lat", "slope")
  # ----------------------------------------------------------------------
  # 1. Trend map
  # ----------------------------------------------------------------------
  sdev  <- sd(df_trend$slope, na.rm = TRUE)
  clamp <- SD_K * sdev
  df_trend$slope_clamped <- pmax(pmin(df_trend$slope, clamp), -clamp)
  p_map <- ggplot(df_trend, aes(lon, lat)) +
    geom_raster(aes(fill = slope_clamped)) +
    geom_sf(
      data = coast,
      color = "black",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    scale_fill_scico(
      palette = "bam",
      name = sprintf("%s trend (%% yr⁻¹, rel. to 1982–2000)", VAR),
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
  global_mean_trend <- with(df_trend,
                            weighted.mean(slope, w = cos(lat * pi / 180), na.rm = TRUE)
  )

  writeLines(
    sprintf(
      "Global mean trend in dropped region: %.6f %% yr-1 (rel. to 1982–2000)",
      global_mean_trend
    ),
    con = file.path(OUTDIR, "summary.txt")
  )
  # ----------------------------------------------------------------------
  # 3. Zonal mean trend
  # ----------------------------------------------------------------------
  df_zonal <- df_trend |>
    mutate(lat_band = round(lat)) |>
    group_by(lat_band) |>
    summarise(
      mean_trend = weighted.mean(slope, cos(lat * pi / 180), na.rm = TRUE)
    )
  p_zonal <- ggplot(df_zonal, aes(mean_trend, lat_band)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_path(linewidth = 0.7, color = "black") +
    scale_y_continuous(labels = lab_deg) +
    labs(
      x = sprintf("%s trend (%% yr⁻¹)", VAR),
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
  r_dropped_ts <- r_geo * dropped_mask
  dropped_ts <- terra::global(r_dropped_ts, "mean", na.rm = TRUE) |>
    as_tibble() |>
    rename(value = 1) |>
    mutate(year = years)

  base_ts <- mean(dropped_ts$value[dropped_ts$year <= 2000], na.rm = TRUE)

  if (!is.finite(base_ts) || abs(base_ts) < EPS) {
    warning("Invalid baseline for dropped-region time series — skipping TS")
  } else {
    dropped_ts <- dropped_ts |> mutate(value = 100 * value / base_ts)
    sig_drop <- trend_test_hac(dropped_ts)

    p_ts <- ggplot(dropped_ts, aes(year, value)) +
      geom_line(linewidth = 0.55) +
      geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.55) +
      labs(
        x = "Year",
        y = sprintf("%s (global mean; %% of 1982–2000)", VAR),
        title = sprintf("%s %s: Annual series (masked-out region)", VAR, METRIC),
        subtitle = sprintf("%s mask | slope = %.3g %% yr⁻¹, p = %.3f",
                           MASK, sig_drop$slope, sig_drop$p)
      ) +
      theme_pub()

    ggsave(file.path(OUTDIR, sprintf("timeseries_dropped_%s_%s.png", METRIC, MASK)),
           p_ts, width = 6.4, height = 4.2, dpi = 330)
  }

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
