## =============================================================================
## 16_compare_masked_unmasked.R
## Compare unmasked (all land) vs masked (natural-only) trends and time series
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(scico)
  library(sf)
  library(rnaturalearth)
  library(scales)
  library(tidyr)
  library(tibble)
  library(here)
})

source(here("R", "utils.R"))
source(here("R", "viz.R"))

BASE_START <- 1982
BASE_END   <- 2000
EPS_BASE   <- 1e-8

compare_georef_masked <- function(VAR,
                                  METRIC,
                                  FILE_GEOREF_SLOPE,
                                  FILE_MASKED_SLOPE,
                                  FILE_GEOREF_DATA,
                                  FILE_MASKED_DATA,
                                  OUTDIR,
                                  SD_K = 2) {

  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

  # Coastline is needed inside plotting functions in THIS scope
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

  # ------------------------------------------------------------------
  # 0) Input checks
  # ------------------------------------------------------------------
  files <- c(FILE_GEOREF_SLOPE, FILE_MASKED_SLOPE, FILE_GEOREF_DATA, FILE_MASKED_DATA)
  if (!all(file.exists(files))) {
    warning("Missing inputs for ", VAR, " / ", METRIC, " in ", OUTDIR, " — skipping.")
    return(invisible(NULL))
  }

  # ------------------------------------------------------------------
  # 1) Relative trend maps (% yr-1) for unmasked vs masked
  # ------------------------------------------------------------------
  r_geo_slope  <- rast(FILE_GEOREF_SLOPE)
  r_mask_slope <- rast(FILE_MASKED_SLOPE)

  r_geo_ts  <- rast(FILE_GEOREF_DATA)
  r_mask_ts <- rast(FILE_MASKED_DATA)

  years_geo  <- BASE_START:(BASE_START + nlyr(r_geo_ts)  - 1L)
  years_mask <- BASE_START:(BASE_START + nlyr(r_mask_ts) - 1L)

  if (!any(years_geo  >= BASE_START & years_geo  <= BASE_END) ||
      !any(years_mask >= BASE_START & years_mask <= BASE_END)) {
    warning("Baseline window not covered for ", VAR, " / ", METRIC, " — skipping.")
    return(invisible(NULL))
  }

  r_geo_base  <- baseline_mean_raster(r_geo_ts,  years_geo,  BASE_START, BASE_END)
  r_mask_base <- baseline_mean_raster(r_mask_ts, years_mask, BASE_START, BASE_END)

  r_geo_rel <- ifel(abs(r_geo_base) < EPS_BASE, NA_real_, 100 * (r_geo_slope  / r_geo_base))
  r_mask_rel <- ifel(abs(r_mask_base) < EPS_BASE, NA_real_, 100 * (r_mask_slope / r_mask_base))

  df_geo  <- as.data.frame(r_geo_rel,  xy = TRUE, na.rm = TRUE)
  df_mask <- as.data.frame(r_mask_rel, xy = TRUE, na.rm = TRUE)
  colnames(df_geo)  <- c("lon", "lat", "slope")
  colnames(df_mask) <- c("lon", "lat", "slope")

  p_map_geo <- plot_map_slope(
    df_geo,
    title_lab = paste0(VAR, " ", METRIC, " trend — all land"),
    SD_K = SD_K
  )
  p_map_mask <- plot_map_slope(
    df_mask,
    title_lab = paste0(VAR, " ", METRIC, " trend — natural-only"),
    SD_K = SD_K
  )

  ggsave(file.path(OUTDIR, sprintf("%s_%s_trend_map_unmasked.png", VAR, METRIC)),
         p_map_geo, width = 7.2, height = 3.8, dpi = 330)
  ggsave(file.path(OUTDIR, sprintf("%s_%s_trend_map_masked.png", VAR, METRIC)),
         p_map_mask, width = 7.2, height = 3.8, dpi = 330)

  # Merge pixelwise
  df <- inner_join(df_geo, df_mask, by = c("lon", "lat"), suffix = c("_geo", "_mask"))

  # ------------------------------------------------------------------
  # 2) Global annual time series (% of baseline) + HAC tests
  # ------------------------------------------------------------------
  glob_geo <- terra::global(r_geo_ts, "mean", na.rm = TRUE) |>
    as_tibble() |>
    rename(value = 1) |>
    mutate(year = years_geo, type = "All land (unmasked)")

  glob_mask <- terra::global(r_mask_ts, "mean", na.rm = TRUE) |>
    as_tibble() |>
    rename(value = 1) |>
    mutate(year = years_mask, type = "Natural-only (masked)")

  base_geo <- mean(glob_geo$value[glob_geo$year >= BASE_START & glob_geo$year <= BASE_END], na.rm = TRUE)
  base_mask <- mean(glob_mask$value[glob_mask$year >= BASE_START & glob_mask$year <= BASE_END], na.rm = TRUE)

  if (!is.finite(base_geo) || abs(base_geo) < EPS_BASE ||
      !is.finite(base_mask) || abs(base_mask) < EPS_BASE) {
    warning("Invalid baseline mean in ", OUTDIR, " — skipping time series.")
  } else {
    glob_geo  <- glob_geo  |> mutate(value = 100 * value / base_geo)
    glob_mask <- glob_mask |> mutate(value = 100 * value / base_mask)

    sig_geo  <- trend_test_hac(glob_geo)
    sig_mask <- trend_test_hac(glob_mask)

    glob_long <- bind_rows(glob_geo, glob_mask)

    p_ts <- ggplot(glob_long, aes(year, value, colour = type)) +
      geom_line(linewidth = 0.6) +
      scale_colour_manual(
        values = c("All land (unmasked)" = "#1b9e77",
                   "Natural-only (masked)" = "#d95f02"),
        name = NULL
      ) +
      labs(
        x = "Year",
        y = sprintf("%s (%s; %% of %d–%d mean)", VAR, METRIC, BASE_START, BASE_END),
        title = sprintf("%s %s: global mean comparison", VAR, METRIC),
        subtitle = sprintf(
          "HAC slopes (%% yr\u207B\u00B9): unmasked %.3g (p=%.3f), masked %.3g (p=%.3f)",
          sig_geo$slope, sig_geo$p, sig_mask$slope, sig_mask$p
        )
      ) +
      theme_pub() +
      theme(legend.position = "bottom")

    ggsave(file.path(OUTDIR, sprintf("%s_%s_global_timeseries_geo_vs_mask.png", VAR, METRIC)),
           p_ts, width = 6.8, height = 4.4, dpi = 330)
  }

  # ------------------------------------------------------------------
  # 3) Scatterplot: pixelwise trends (masked vs unmasked)
  # ------------------------------------------------------------------
  df_stats <- df |>
    filter(is.finite(slope_geo), is.finite(slope_mask))

  if (nrow(df_stats) > 10) {
    lmfit <- lm(slope_mask ~ slope_geo, data = df_stats)
    r_pearson <- cor(df_stats$slope_geo, df_stats$slope_mask, use = "complete.obs")
    slope_lm  <- unname(coef(lmfit)[2])
    r2_lm     <- summary(lmfit)$r.squared
    N         <- nrow(df_stats)
    sign_agree <- mean(sign(df_stats$slope_geo) == sign(df_stats$slope_mask)) * 100

    stats_label <- sprintf(
      "r = %.2f | R\u00B2 = %.2f | slope = %.2f | sign agree = %.1f%% | N = %d",
      r_pearson, r2_lm, slope_lm, sign_agree, N
    )

    p_sc <- ggplot(df_stats, aes(slope_geo, slope_mask)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = "grey50", linewidth = 0.4) +
      geom_point(alpha = 0.25, size = 0.3, colour = "black") +
      geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.7) +
      labs(
        x = "All land (unmasked) trend (% yr⁻¹)",
        y = "Natural-only (masked) trend (% yr⁻¹)",
        title = sprintf("%s %s: pixel-wise trend comparison", VAR, METRIC),
        subtitle = stats_label
      ) +
      theme_pub()

    ggsave(file.path(OUTDIR, sprintf("%s_%s_scatter_geo_vs_mask.png", VAR, METRIC)),
           p_sc, width = 6.2, height = 5.6, dpi = 330)
  }

  # ------------------------------------------------------------------
  # 4) Zonal means: trend vs latitude (area-weighted)
  # ------------------------------------------------------------------
  df_zon <- df |>
    mutate(
      lat_band = floor(lat),
      w = cos(lat * pi / 180)
    ) |>
    group_by(lat_band) |>
    summarise(
      geo  = weighted.mean(slope_geo,  w, na.rm = TRUE),
      mask = weighted.mean(slope_mask, w, na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_longer(cols = c("geo", "mask"), names_to = "type", values_to = "trend") |>
    mutate(type = recode(type,
                         geo = "All land (unmasked)",
                         mask = "Natural-only (masked)"))

  p_zonal <- ggplot(df_zon, aes(trend, lat_band, colour = type)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_path(linewidth = 0.8) +
    scale_y_continuous(
      limits = c(-90, 90),
      breaks = seq(-90, 90, 30),
      labels = lab_deg
    ) +
    scale_colour_manual(
      values = c("All land (unmasked)" = "#1b9e77",
                 "Natural-only (masked)" = "#d95f02"),
      name = NULL
    ) +
    labs(
      x = "Trend (% yr⁻¹)",
      y = "Latitude (°)",
      title = sprintf("%s %s: zonal trends", VAR, METRIC),
      subtitle = "Area-weighted 1° zonal means"
    ) +
    theme_pub() +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTDIR, sprintf("%s_%s_zonal_geo_vs_mask.png", VAR, METRIC)),
         p_zonal, width = 6.4, height = 5.4, dpi = 330)

  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

VARS    <- c("LAI", "FPAR")
METRICS <- c("yearmean", "yearmax")
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
