## =============================================================================
## 15_analyse_masked_trends.R
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(scico)
  library(sf)
})

source(here("R", "utils.R"))
source(here("R", "viz.R"))

BASE_START <- 1982
BASE_END   <- 2000
MIN_BASE_YEARS <- 10L

analyse_masked_trends <- function(ROOT_DIR, VAR, OUTDIR, SD_K = 2) {

  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

  metrics <- c("yearmean", "yearmax")

  for (M in metrics) {

    message("Processing ", VAR, " / ", M, " — ", basename(ROOT_DIR))

    f_data <- file.path(ROOT_DIR, sprintf("%s_%s_0p25.nc", VAR, M))
    f_year <- file.path(ROOT_DIR,
                        sprintf("%s_%s_trend_slope_peryear_0p25.nc", VAR, M))

    if (!file.exists(f_data) || !file.exists(f_year)) {
      warning("Missing input for ", VAR, " / ", M, ", skipping.")
      next
    }

    # ------------------------------------------------------------------
    # Load annual series
    # ------------------------------------------------------------------
    r_series <- rast(f_data)
    years    <- BASE_START:(BASE_START + nlyr(r_series) - 1)

    if (!any(years >= BASE_START & years <= BASE_END)) {
      warning("Baseline window not covered for ", basename(f_data))
      next
    }

    # ------------------------------------------------------------------
    # Global time series (area-mean)
    # ------------------------------------------------------------------
    w <- cos(yFromRow(r_series, 1:nrow(r_series)) * pi/180)
    w_r <- init(r_series[[1]], "y")
    values(w_r) <- rep(w, each = ncol(r_series))
    glob <- terra::global(r_series, "mean", w = w_r, na.rm = TRUE) |>
      as_tibble() |>
      rename(value = 1) |>
      mutate(year = years)

    base_mean <- mean(
      glob$value[glob$year >= BASE_START & glob$year <= BASE_END],
      na.rm = TRUE
    )

    if (!is.finite(base_mean) || abs(base_mean) < 1e-8) {
      warning("Invalid baseline mean for ", basename(ROOT_DIR), " / ", M)
      next
    }

    glob <- glob |> mutate(value = 100 * value / base_mean)
    sig <- trend_test_hac(glob)
    p_ts <- ggplot(glob, aes(year, value)) +
      geom_line(linewidth = 0.55, color = "black") +
      geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.55) +
      labs(
        x = "Year",
        y = sprintf("%s (%% of %d–%d mean)", VAR, BASE_START, BASE_END),
        title = sprintf(
          "%s: Global mean over natural vegetation — %s",
          VAR, M
        ),
        subtitle = sprintf(
          "%s | slope = %.3g %% yr⁻¹, p = %.3f (HAC)",
          basename(ROOT_DIR), sig$slope, sig$p
        )
      ) +
      theme_pub()

    ggsave(
      file.path(OUTDIR, sprintf("%s_%s_timeseries.png", VAR, M)),
      p_ts,
      width = 6.5,
      height = 4.2,
      dpi = 330
    )

    # ------------------------------------------------------------------
    # Pixel-wise relative trends
    # ------------------------------------------------------------------
    r_slope <- rast(f_year)

    idx_base <- which(years >= BASE_START & years <= BASE_END)
    r_base   <- mean(r_series[[idx_base]], na.rm = TRUE)

    n_ok_base <- terra::app(
      r_series[[idx_base]],
      fun = function(v) sum(is.finite(v))
    )

    r_rel <- ifel(
      (abs(r_base) < 1e-8) | (n_ok_base < MIN_BASE_YEARS),
      NA_real_,
      100 * (r_slope / r_base)
    )
    stopifnot(isTRUE(all.equal(ext(r_rel), ext(r_slope))))

    df_year <- as.data.frame(r_rel, xy = TRUE, na.rm = TRUE)
    colnames(df_year) <- c("lon", "lat", "slope")

    # ------------------------------------------------------------------
    # Trend map
    # ------------------------------------------------------------------
    p_map <- plot_map_slope(df_year, VAR, SD_K = SD_K)

    ggsave(
      file.path(OUTDIR, sprintf("%s_%s_map_year.png", VAR, M)),
      p_map,
      width = 7.2,
      height = 3.8,
      dpi = 330
    )

    # ------------------------------------------------------------------
    # Zonal mean profile
    # ------------------------------------------------------------------
    df_zonal <- df_year |>
      mutate(
        lat_band = floor(lat),
        w = cos(lat * pi / 180)
      ) |>
      group_by(lat_band) |>
      summarise(
        slope_mean = weighted.mean(slope, w, na.rm = TRUE),
        .groups = "drop"
      )

    p_zonal <- ggplot(df_zonal, aes(slope_mean, lat_band)) +
      geom_vline(xintercept = 0, color = "grey60", linewidth = 0.35) +
      geom_path(color = "black", linewidth = 0.55) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
      scale_y_continuous(
        limits = c(-90, 90),
        breaks = seq(-90, 90, 30),
        labels = lab_deg
      ) +
      labs(
        x = sprintf("%s trend (%% yr⁻¹)", VAR),
        y = "Latitude (°)",
        title = sprintf("%s: Zonal trend — %s", VAR, M),
        subtitle = basename(ROOT_DIR)
      ) +
      theme_pub()

    ggsave(
      file.path(OUTDIR, sprintf("%s_%s_zonal_year.png", VAR, M)),
      p_zonal,
      width = 5.3,
      height = 4.8,
      dpi = 330
    )
  }

  message("Completed: ", basename(ROOT_DIR))
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------
TAU <- c("tau_0.05", "tau_0.1", "tau_0.2")

for (tau in TAU) {

  BASE <- here(sprintf("output/%s/eval", tau))

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
