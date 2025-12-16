## =============================================================================
## 13_analyse_all_georef_0p25.R — Scientific diagnostics for LAI/FPAR (0.25°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(rnaturalearth)
  library(sf)
  library(here)
  library(scico)
  library(scales)
  library(lmtest)
  library(sandwich)
})

source(here("R", "utils.R"))
source(here("R", "viz.R"))

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# -----------------------------------------------------------------------------
VAR <- "LAI"   # LAI or "FPAR"
ROOT <- here::here()
SD_K <- 2      # clamp at ±K standard deviations
BASE_START <- 1982
BASE_END   <- 2000


metrics <- c("yearmean", "yearmax")

input_dir <- file.path(ROOT, "analysis", "unmasked", "0p25")
fig_dir   <- file.path(ROOT, "analysis", "trends_unmasked", VAR)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
lab_deg <- label_deg()


# -----------------------------------------------------------------------------
# Helper
# -----------------------------------------------------------------------------

plot_trend_map <- function(df, metric, trend_label, var_name, SD_K = 2) {
  sdev  <- sd(df$slope, na.rm = TRUE)
  clamp <- SD_K * sdev
  df$slope_clamped <- pmax(pmin(df$slope, clamp), -clamp)

  ggplot(df, aes(lon, lat)) +
    geom_raster(aes(fill = slope_clamped)) +
    geom_sf(data = coast, colour = "black", linewidth = 0.15, inherit.aes = FALSE) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 60), labels = lab_deg) +
    scale_y_continuous(breaks = seq(-90, 90, 30), labels = lab_deg) +
    scale_fill_scico(
      palette = "bam",
      name = paste0(var_name, " trend (% yr⁻¹, baseline-normalised)"),
      limits = c(-clamp, clamp),
      oob    = scales::squish
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = sprintf("%s: %s trend (%s)", var_name, metric, trend_label)
    ) +
    theme_pub()
}


plot_zonal <- function(df_zonal, metric, trend_label, var_name) {
  ggplot(df_zonal, aes(slope_mean, lat_band)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.35) +
    geom_path(linewidth = 0.7, color = "black") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "red") +
    scale_y_continuous(
      limits = c(-90, 90),
      breaks = seq(-90, 90, 30),
      labels = lab_deg
    ) +
    labs(
      x = sprintf("%s trend (%s)", var_name, trend_label),
      y = "Latitude (°)",
      title = sprintf("%s: zonal mean %s trend", var_name, metric),
      subtitle = "1° latitude bands"
    ) +
    theme_pub()
}

# -----------------------------------------------------------------------------
# Loop over metrics
# -----------------------------------------------------------------------------
for (M in metrics) {
  cat("\n============================================================\n")
  cat("Processing metric:", M, "\n")
  cat("============================================================\n")

  # --- File discovery ---------------------------------------------------------
  file_data   <- file.path(input_dir, sprintf("%s_georef_%s_0p25.nc", VAR, M))
  file_slope_year <- file.path(input_dir,
                               sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", VAR, M))

  if (!file.exists(file_data)) {
    warning("Missing data file for metric ", M, ", skipping.")
    next
  }
  if (!file.exists(file_slope_year)) {
    warning("Missing trend files for metric ", M, ", skipping.")
    next
  }

  # --- Load annual series -----------------------------------------------------
  r_series <- rast(file_data)
  years    <- BASE_START:(BASE_START + nlyr(r_series) - 1)
  # --- sanity: baseline years must exist in this file ---
  if (!any(years >= BASE_START & years <= BASE_END)) {
    warning("Baseline window not covered by ", basename(file_data), " — skipping.")
    next
  }

  # --- Global time series -----------------------------------------------------
  glob <- terra::global(r_series, "mean", na.rm = TRUE) |>
    as_tibble() |>
    rename(value = 1) |>
    mutate(year = years)

  # --- relative change (% of baseline mean) ---
  base_mean <- mean(glob$value[glob$year >= BASE_START & glob$year <= BASE_END], na.rm = TRUE)
  glob <- glob |> mutate(value = 100 * value / base_mean)
  sig  <- trend_test_hac(glob)

  p_ts <- ggplot(glob, aes(year, value)) +
    geom_line(linewidth = 0.55, color = "black") +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "red",
      linewidth = 0.55
    ) +
    labs(
      x = "Year",
      y = sprintf("%s (%% of %d–%d mean)", VAR, BASE_START, BASE_END),
      title = paste0(VAR, ": Global annual mean — ", M),
      subtitle = sprintf(
        "Relative slope = %.3g %% yr⁻¹, p = %.3f (HAC)",
        sig$slope, sig$p
      )
    ) +
    theme_pub()

  ggsave(
    file.path(fig_dir, sprintf("%s_%s_global_timeseries.png", VAR, M)),
    p_ts,
    width = 6.5,
    height = 4.2,
    dpi = 330
  )

  # --- Load slope data --------------------------------------------------------
  r_year <- rast(file_slope_year)
  # Require at least N valid baseline years per pixel to define % trends robustly
  MIN_BASE_YEARS <- 10L
  idx_base <- which(years >= BASE_START & years <= BASE_END)
  n_ok_base <- terra::app(r_series[[idx_base]], fun = function(v) sum(is.finite(v)))

  # --- relative trend (% per year) ---
  r_base <- baseline_mean_raster(r_series, years, BASE_START, BASE_END)

  # Convert absolute slope (units yr^-1) to relative slope (% yr^-1)
  r_year <- 100 * (r_year / r_base)

  r_year <- ifel(
    (abs(r_base) < 1e-8) | (n_ok_base < MIN_BASE_YEARS),
    NA, r_year
  )


  df_year <- as.data.frame(r_year, xy = TRUE, na.rm = TRUE)
  colnames(df_year) <- c("lon", "lat", "slope")

  # ===================
  # Trend map: per YEAR
  # ===================
  p_map_year <- plot_trend_map(df_year, M, "per year", VAR, SD_K)
  ggsave(
    file.path(fig_dir, sprintf("%s_%s_trend_peryear_map.png", VAR, M)),
    p_map_year,
    width = 7.2,
    height = 3.8,
    dpi = 330
  )

  df_year_zonal <- df_year |>
    mutate(lat_band = floor(lat)) |>
    group_by(lat_band) |>
    summarise(
      slope_mean = weighted.mean(slope, w = cos(lat_band * pi / 180), na.rm = TRUE),
      .groups = "drop"
    )

  p_z_year <- plot_zonal(df_year_zonal, M, "per year", VAR)
  ggsave(
    file.path(fig_dir, sprintf("%s_%s_zonal_peryear.png", VAR, M)),
    p_z_year,
    width = 5.3,
    height = 4.8,
    dpi = 330
  )
}

cat("\nAll analyses completed.\nFigures saved to:\n  ", fig_dir, "\n")

