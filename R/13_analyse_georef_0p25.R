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
  library(ggspatial)
  library(scales)
})

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# -----------------------------------------------------------------------------
VAR <- "LAI"   # LAI or "FPAR"
ROOT <- here::here()
SD_K <- 2      # clamp at ±K standard deviations

metrics <- c("yearmean", "yearmax", "yearamp")

input_dir <- file.path(ROOT, "analysis", "unmasked", "0p25")
fig_dir   <- file.path(ROOT, "analysis", "trends_unmasked", VAR)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
land  <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
grat <- sf::st_graticule(lon = seq(-180, 180, by = 30),
                         lat = seq(-90, 90, by = 30))

# -----------------------------------------------------------------------------
# Helper
# -----------------------------------------------------------------------------
theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey87", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 13, face = "bold"),
      plot.subtitle    = element_text(size = 10),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 9)
    )
}

plot_trend_map <- function(df, metric, trend_label, var_name, SD_K = 2) {
  sdev  <- sd(df$slope, na.rm = TRUE)
  clamp <- SD_K * sdev
  df$slope_clamped <- pmax(pmin(df$slope, clamp), -clamp)
  lab_deg <- scales::label_number(suffix = "°")

  ggplot(df, aes(lon, lat)) +
    # raster field
    geom_raster(aes(fill = slope_clamped)) +
    # coastline overlay
    geom_sf(
      data = coast,
      colour = "black",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    # coordinate system + limits
    coord_sf(
      xlim = range(df$lon, na.rm = TRUE),
      ylim = range(df$lat, na.rm = TRUE),
      expand = FALSE
    ) +
    scale_x_continuous(breaks = seq(-180, 180, 60), labels = lab_deg) +
    scale_y_continuous(breaks = seq(-90, 90, 30), labels = lab_deg) +
    scale_fill_scico(
      palette = "bam",
      name   = paste0(var_name, " trend (", trend_label, ")"),
      limits = c(-clamp, clamp),
      oob    = scales::squish
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title    = paste0(var_name, ": ", metric, " trend (", trend_label, ")"),
    )
}

plot_zonal <- function(df_zonal, metric, trend_label, var_name) {
  ggplot(df_zonal, aes(slope_mean, lat_band)) +
    geom_vline(xintercept = 0,
               color = "grey60",
               linewidth = 0.35) +
    geom_path(color = "black", linewidth = 0.55) +
    labs(
      x = paste0(var_name, " trend (", trend_label, ")"),
      y = "Latitude (°)",
      title = paste0(var_name, ": Zonal mean trend — ", metric),
      subtitle = "Mean trend for 1° latitude bands"
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
  years    <- 1982:(1982 + nlyr(r_series) - 1)

  # --- Global time series -----------------------------------------------------
  glob <- terra::global(r_series, "mean", na.rm = TRUE) |>
    as_tibble() |>
    rename(value = 1) |>
    mutate(year = years)

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
      y = paste0(VAR, " (annual global mean)"),
      title = paste0(VAR, ": Global annual mean — ", M),
      subtitle = sprintf("Period: %d–%d (unmasked 0.25°)", years[1], years[length(years)])
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
  ref <- rast(file.path(ROOT, "src/ref_0p25.nc"))
  r_year <- terra::project(r_year, ref)
  crs(r_year) <- "EPSG:4326"
  ext(r_year) <- ext(-180, 180, -90, 90)
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
    summarise(slope_mean = mean(slope, na.rm = TRUE),
              .groups = "drop")
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

