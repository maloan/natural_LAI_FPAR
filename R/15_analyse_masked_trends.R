## =============================================================================
## 15_analyse_masked_trends.R
## =============================================================================

analyse_masked_trends <- function(ROOT_DIR, VAR, OUTDIR, SD_K = 2) {
  suppressPackageStartupMessages({
    library(terra)
    library(dplyr)
    library(ggplot2)
    library(scico)
    library(sf)
    library(rnaturalearth)
    library(scales)
  })
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
  theme_pub <- function() {
    theme_bw(base_size = 12) +
      theme(
        panel.grid.major = element_line(color = "grey87", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.text  = element_text(size = 9)
      )
  }
  lab_deg <- scales::label_number(suffix = "°")
  plot_trend_map <- function(df, metric, trend_label, var_name, SD_K = 2) {
    sdev  <- sd(df$slope, na.rm = TRUE)
    clamp <- SD_K * sdev
    df$slope_clamped <- pmax(pmin(df$slope, clamp), -clamp)

    ggplot(df, aes(lon, lat)) +
      geom_raster(aes(fill = slope_clamped)) +
      geom_sf(
        data = coast,
        colour = "black",
        linewidth = 0.2,
        inherit.aes = FALSE
      ) +
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
        title    = paste0(var_name, ": ", metric, " trend (", trend_label, ")")
      ) +
      theme_pub()
  }

  metrics <- c("yearmean", "yearmax")
  for (M in metrics) {
    f_data <- file.path(ROOT_DIR, sprintf("%s_%s_0p25.nc", VAR, M))
    f_year <- file.path(ROOT_DIR,
                        sprintf("%s_%s_trend_slope_peryear_0p25.nc", VAR, M))
    # -------------------------------------------------------------
    # Global time series
    # -------------------------------------------------------------
    r_series <- rast(f_data)
    years <- seq(1982, 1982 + nlyr(r_series) - 1)
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
        y = paste0(VAR, " (global mean)"),
        title = paste0(VAR, ": Global annual mean — ", M),
        subtitle = sprintf("%s", basename(ROOT_DIR))
      ) +
      theme_pub()

    ggsave(
      file.path(OUTDIR, sprintf("%s_%s_timeseries.png", VAR, M)),
      p_ts,
      width = 6.5,
      height = 4.2,
      dpi = 330
    )

    # -------------------------------------------------------------
    # Load slope rasters and convert to df
    # -------------------------------------------------------------
    r_year <- rast(f_year)
    df_year <- as.data.frame(r_year, xy = TRUE, na.rm = TRUE)
    colnames(df_year) <- c("lon", "lat", "slope")

    # -------------------------------------------------------------
    # Trend maps
    # -------------------------------------------------------------
    p_year <- plot_trend_map(df_year, M, "per year", VAR, SD_K)
    ggsave(
      file.path(OUTDIR, sprintf("%s_%s_map_year.png", VAR, M)),
      p_year,
      width = 7.2,
      height = 3.8,
      dpi = 330
    )

    # -------------------------------------------------------------
    # Zonal mean profiles
    # -------------------------------------------------------------
    df_year_z <- df_year |>
      mutate(lat_band = floor(lat)) |>
      group_by(lat_band) |>
      summarise(slope_mean = mean(slope, na.rm = TRUE))

    p_z_year <- ggplot(df_year_z, aes(slope_mean, lat_band)) +
      geom_vline(
        xintercept = 0,
        color = "grey60",
        linewidth = 0.35
      ) +
      geom_path(color = "black", linewidth = 0.55) +
      labs(
        x = paste0(VAR, " trend / year"),
        y = "Latitude (°)",
        title = paste0(VAR, ": Zonal trend — ", M),
        subtitle = basename(ROOT_DIR)
      ) +
      theme_pub()
    ggsave(
      file.path(OUTDIR, sprintf("%s_%s_zonal_year.png", VAR, M)),
      p_z_year,
      width = 5.3,
      height = 4.8,
      dpi = 330
    )
  }

  message("\nCompleted: ", basename(ROOT_DIR))
}
# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------
TAU <- c("tau_0.05", "tau_0.1", "tau_0.2")

for (tau in TAU) {
  BASE <- here::here(sprintf("output/%s/eval", tau))
  analyse_masked_trends(
    ROOT_DIR = file.path(BASE, "trend_FPAR_CCI"),
    VAR = "FPAR",
    OUTDIR = sprintf("analysis/trends_masked/%s/FPAR_CCI", tau)
  )

  analyse_masked_trends(
    ROOT_DIR = file.path(BASE, "trend_FPAR_GLC"),
    VAR = "FPAR",
    OUTDIR = sprintf("analysis/trends_masked/%s/FPAR_GLC", tau)
  )

  analyse_masked_trends(
    ROOT_DIR = file.path(BASE, "trend_LAI_CCI"),
    VAR = "LAI",
    OUTDIR = sprintf("analysis/trends_masked/%s/LAI_CCI", tau)
  )

  analyse_masked_trends(
    ROOT_DIR = file.path(BASE, "trend_LAI_GLC"),
    VAR = "LAI",
    OUTDIR = sprintf("analysis/trends_masked/%s/LAI_GLC", tau)
  )
}
