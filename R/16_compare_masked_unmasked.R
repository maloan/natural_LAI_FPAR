
## =============================================================================
## 16_compare_masked_unmasked.R
## =============================================================================
compare_georef_masked <- function(VAR,
                                  METRIC,
                                  FILE_GEOREF_SLOPE,
                                  FILE_MASKED_SLOPE,
                                  FILE_GEOREF_DATA,
                                  FILE_MASKED_DATA,
                                  OUTDIR,
                                  SD_K = 2) {
  suppressPackageStartupMessages({
    library(terra)
    library(dplyr)
    library(ggplot2)
    library(scico)
    library(sf)
    library(rnaturalearth)
    library(scales)
    library(tidyr)
  })

  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
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

  lab_deg <- scales::label_number(suffix = "°")
  plot_trend_map <- function(df, title_lab, SD_K = 2) {
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
      scale_fill_scico(palette = "bam",
                       name    = title_lab,
                       oob     = scales::squish) +
      labs(x = "Longitude", y = "Latitude") +
      theme_pub()
  }

  # ------------------------------------------------------------------
  # 1) Load slope rasters (per year trends) + trend maps
  # ------------------------------------------------------------------
  r_geo  <- rast(FILE_GEOREF_SLOPE)
  r_mask <- rast(FILE_MASKED_SLOPE)
  df_geo  <- as.data.frame(r_geo, xy = TRUE, na.rm = TRUE)
  df_mask <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  colnames(df_geo)  <- c("lon", "lat", "slope")
  colnames(df_mask) <- c("lon", "lat", "slope")
  # maps: unmasked vs masked
  p_map_geo <- plot_trend_map(df_geo,
                              paste0(VAR, " ", METRIC, " trend (per year, all land)"),
                              SD_K)
  p_map_mask <- plot_trend_map(df_mask,
                               paste0(VAR, " ", METRIC, " trend (per year, natural-only)"),
                               SD_K)
  ggsave(
    file.path(
      OUTDIR,
      sprintf("%s_%s_trend_map_unmasked.png", VAR, METRIC)
    ),
    p_map_geo,
    width = 7.2,
    height = 3.8,
    dpi = 330
  )
  ggsave(
    file.path(OUTDIR, sprintf(
      "%s_%s_trend_map_masked.png", VAR, METRIC
    )),
    p_map_mask,
    width = 7.2,
    height = 3.8,
    dpi = 330
  )
  # Merge (slope_geo, slope_mask)
  df <- dplyr::inner_join(df_geo,
                          df_mask,
                          by = c("lon", "lat"),
                          suffix = c("_geo", "_mask"))

  # ------------------------------------------------------------------
  # 2) Global annual time series: all land vs natural-only
  # ------------------------------------------------------------------
  r_geo_ts  <- rast(FILE_GEOREF_DATA)
  r_mask_ts <- rast(FILE_MASKED_DATA)
  years_geo  <- 1982:(1982 + nlyr(r_geo_ts)  - 1L)
  years_mask <- 1982:(1982 + nlyr(r_mask_ts) - 1L)
  glob_geo <- terra::global(r_geo_ts, "mean", na.rm = TRUE) |>
    as_tibble() |>
    dplyr::rename(value = 1) |>
    dplyr::mutate(year = years_geo, type = "All land (unmasked)")
  glob_mask <- terra::global(r_mask_ts, "mean", na.rm = TRUE) |>
    as_tibble() |>
    dplyr::rename(value = 1) |>
    dplyr::mutate(year = years_mask, type = "Natural-only (masked)")
  glob_long <- dplyr::bind_rows(glob_geo, glob_mask)
  p_ts <- ggplot(glob_long, aes(year, value, colour = type)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(
      values = c(
        "All land (unmasked)"   = "#1b9e77",
        "Natural-only (masked)" = "#d95f02"
      ),
      name = NULL
    ) +
    labs(
      x = "Year",
      y = paste0(VAR, " (annual global ", METRIC, ")"),
      title = paste0(VAR, " ", METRIC, ": global annual mean"),
      subtitle = "Comparison of all land (unmasked) vs natural-only (masked)"
    ) +
    theme_pub() +
    theme(legend.position = "bottom",
          legend.margin   = margin(t = 2, b = 2))
  ggsave(
    file.path(
      OUTDIR,
      sprintf("%s_%s_global_timeseries_geo_vs_mask.png", VAR, METRIC)
    ),
    p_ts,
    width = 6.8,
    height = 4.4,
    dpi = 330
  )
  # ------------------------------------------------------------------
  # 3) Scatterplot: pixelwise slopes (masked vs unmasked)
  # ------------------------------------------------------------------
  df_stats <- df |>
    dplyr::filter(is.finite(slope_geo), is.finite(slope_mask))
  lmfit <- lm(slope_mask ~ slope_geo, data = df_stats)
  r_pearson <- stats::cor(df_stats$slope_geo, df_stats$slope_mask)
  slope_lm  <- coef(lmfit)[2]
  r2_lm     <- summary(lmfit)$r.squared
  N         <- nrow(df_stats)
  sign_agree <- mean(sign(df_stats$slope_geo) == sign(df_stats$slope_mask)) * 100
  stats_label <- sprintf(
    "r = %.2f | R² = %.2f | slope = %.2f\nSign agreement: %.1f%% | N = %d",
    r_pearson,
    r2_lm,
    slope_lm,
    sign_agree,
    N
  )
  sdev_sc  <- sd(df$slope, na.rm = TRUE)
  # clamp_sc <- SD_K * sdev_sc
  # lim_sc   <- c(-clamp_sc, clamp_sc)
  p_sc <- ggplot(df_stats, aes(slope_geo, slope_mask)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.4
    ) +
    geom_point(alpha = 0.25,
               size = 0.3,
               colour = "black") +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "red",
      linewidth = 0.7
    ) +
    # coord_equal(xlim = lim_sc,
    #             ylim = lim_sc,
    #             expand = TRUE) +
    labs(
      x = "All land (unmasked) trend (per year)",
      y = "Natural-only (masked) trend (per year)",
      title = paste0(VAR, " ", METRIC, ": pixel-wise trends"),
      subtitle = paste0(
        "Unmasked = includes cropland, urban, pasture; ",
        "Masked = only natural vegetation\n",
        stats_label
      )
    ) +
    theme_pub()
  ggsave(
    file.path(
      OUTDIR,
      sprintf("%s_%s_scatter_geo_vs_mask.png", VAR, METRIC)
    ),
    p_sc,
    width = 6.2,
    height = 5.6,
    dpi = 330
  )

  # ------------------------------------------------------------------
  # 4) Zonal means: trend vs latitude (lat on y axis)
  # ------------------------------------------------------------------
  df$lat_band <- floor(df$lat)
  zon <- df |>
    dplyr::group_by(lat_band) |>
    dplyr::summarise(
      geo  = mean(slope_geo, na.rm = TRUE),
      mask = mean(slope_mask, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = c("geo", "mask"),
      names_to = "type",
      values_to = "trend"
    ) |>
    dplyr::mutate(type = dplyr::recode(type, geo  = "All land (unmasked)", mask = "Natural-only (masked)"))

  p_zonal <- ggplot(zon, aes(trend, lat_band, colour = type)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_path(linewidth = 0.8) +
    scale_y_continuous(labels = lab_deg) +
    scale_colour_manual(
      values = c(
        "All land (unmasked)"   = "#1b9e77",
        "Natural-only (masked)" = "#d95f02"
      ),
      name = NULL
    ) +
    labs(
      x = "Trend (per year)",
      y = "Latitude (°)",
      title = paste0(VAR, " ", METRIC, ": zonal trends"),
      subtitle = "Comparison of all land vs natural-only (masked) 1° zonal means"
    ) +
    theme_pub() +
    theme(legend.position = "bottom",
          legend.margin   = margin(t = 2, b = 2))

  ggsave(
    file.path(OUTDIR, sprintf(
      "%s_%s_zonal_geo_vs_mask.png", VAR, METRIC
    )),
    p_zonal,
    width = 6.4,
    height = 5.4,
    dpi = 330
  )
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

VARS    <- c("LAI", "FPAR")
METRICS <- c("yearmean", "yearmax", "yearamp")
MASKS   <- c("CCI", "GLC")
TAU     <- c("tau_0.05", "tau_0.1", "tau_0.2")

for (tau in TAU) {
  for (var in VARS) {
    for (met in METRICS) {
      for (mask in MASKS) {
        compare_georef_masked(
          VAR    = var,
          METRIC = met,
          FILE_GEOREF_SLOPE = file.path(
            "analysis/unmasked/0p25",
            sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, met)
          ),
          FILE_MASKED_SLOPE = file.path(
            sprintf("output/%s/eval", tau),
            sprintf("trend_%s_%s", var, mask),
            sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, met)
          ),
          FILE_GEOREF_DATA = file.path(
            "analysis/unmasked/0p25",
            sprintf("%s_georef_%s_0p25.nc", var, met)
          ),
          FILE_MASKED_DATA = file.path(
            sprintf("output/%s/eval", tau),
            sprintf("trend_%s_%s", var, mask),
            sprintf("%s_%s_0p25.nc", var, met)
          ),
          OUTDIR = file.path(
            sprintf("analysis/comparison_masked_unmasked/%s", tau),
            sprintf("%s_%s", var, mask)
          )
        )
      }
    }
  }
}
