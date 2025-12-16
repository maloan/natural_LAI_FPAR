
## =============================================================================
## 19_global_timeseries_masks.R — Global annual mean LAI/FPAR (masked vs unmasked)
## =============================================================================
suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(readr)
})

ROOT    <- here::here()
IN_GEO  <- file.path(ROOT, "analysis", "unmasked", "0p25")
OUTDIR  <- file.path(ROOT, "analysis", "global_timeseries_masks")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

VARS    <- c("LAI", "FPAR")
METRICS <- c("yearmean", "yearmax")
MASKS   <- c("CCI", "GLC")
TAU     <- c("tau_0.05", "tau_0.1", "tau_0.2")
BASE_START <- 1982
BASE_END   <- 2000

# ----------------------------------------------------------------------
# Helper
# ----------------------------------------------------------------------
theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey87", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 13, face = "bold"),
      plot.subtitle    = element_text(size = 10),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 9),
      legend.position  = "bottom"
    )
}
# ------------------------------------------------------------------------------
# Trend helpers
# ------------------------------------------------------------------------------

trend_test_hac <- function(df, y = "value", x = "year") {
  suppressPackageStartupMessages({
    library(lmtest)
    library(sandwich)
  })
  fit <- lm(reformulate(x, y), data = df)
  ct  <- lmtest::coeftest(fit, vcov. = sandwich::NeweyWest(fit))
  tibble::tibble(
    slope = unname(coef(fit)[2]),
    p     = ct[2, 4],
    r2    = summary(fit)$r.squared
  )
}


# ----------------------------------------------------------------------
# Main loop over VAR × METRIC
# ----------------------------------------------------------------------
for (var in VARS) {
  message("==== ", var, " ====")

  ts_list <- list()

  for (met in METRICS) {

    # -------------------------------
    # 1) Unmasked global time series
    # -------------------------------
    f_geo <- file.path(IN_GEO, sprintf("%s_georef_%s_0p25.nc", var, met))
    if (file.exists(f_geo)) {
      r_geo <- rast(f_geo)
      years <- 1982:(1982 + nlyr(r_geo) - 1L)

      df_geo <- terra::global(r_geo, "mean", na.rm = TRUE) |>
        as_tibble() |>
        rename(value = 1) |>
        mutate(
          year   = years,
          metric = met,
          domain = "Unmasked (all land)"
        )

      base <- mean(df_geo$value[df_geo$year >= BASE_START & df_geo$year <= BASE_END], na.rm = TRUE)
      if (!is.finite(base) || abs(base) < 1e-8) next
      df_geo <- df_geo |> mutate(value = 100 * value / base)

      ts_list[[length(ts_list) + 1]] <- df_geo
    }

    # -------------------------------
    # 2) Masked global time series
    # -------------------------------
    for (tau in TAU) {
      for (mask in MASKS) {

        f_mask <- file.path(
          ROOT,
          "output",
          tau,
          "eval",
          sprintf("trend_%s_%s", var, mask),
          sprintf("%s_%s_0p25.nc", var, met)
        )
        if (!file.exists(f_mask)) next

        r_mask <- rast(f_mask)
        years  <- 1982:(1982 + nlyr(r_mask) - 1L)

        df_mask <- terra::global(r_mask, "mean", na.rm = TRUE) |>
          as_tibble() |>
          rename(value = 1) |>
          mutate(
            year   = years,
            metric = met,
            domain = sprintf("%s (τ=%s)", mask, sub("tau_", "", tau))
          )

        base_mask <- mean(
          df_mask$value[df_mask$year >= BASE_START & df_mask$year <= BASE_END],
          na.rm = TRUE
        )
        if (!is.finite(base_mask) || abs(base_mask) < 1e-8) next
        df_mask <- df_mask |> mutate(value = 100 * value / base_mask)

        ts_list[[length(ts_list) + 1]] <- df_mask
      }
    }
  }


    df_ts <- bind_rows(ts_list)
    df_ts <- df_ts |> mutate(domain = factor(domain, levels = unique(domain)))
    sig_table <- df_ts |>
      group_by(domain, metric) |>
      group_modify(~ trend_test_hac(.x)) |>
      ungroup() |>
      mutate(VAR = var)

    write_csv(
      sig_table,
      file.path(OUTDIR, sprintf("%s_global_trend_significance.csv", var))
    )


    # ------------------------------------------------------------------
    # Plot: all domains overlapped
    # ------------------------------------------------------------------
    p <- ggplot(df_ts, aes(year, value, colour = domain)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~ metric, ncol = 1, scales = "free_y") +
      labs(
        x = "Year",
        y = sprintf("%s (%% of 1982–2000 mean)", var),
        title = sprintf("%s: global annual trends (masked vs unmasked)", var),
        subtitle = "Each panel is a metric; HAC-tested trends in the table"
      ) +
      theme_pub()


    outf <- file.path(
      OUTDIR,
      sprintf("%s_global_timeseries_yearmean_yearmax_all_masks.png", var)
    )

    ggsave(outf, p, width = 7.8, height = 5.0, dpi = 330)

    # Also save the table for later analysis
    write_csv(df_ts, file.path(
      OUTDIR,
      sprintf("%s_global_timeseries_all_masks.csv", var)
    ))

    message("Saved: ", outf)
}

cat("Finished global time series (masked vs unmasked).\nOutput in:\n  ",
    OUTDIR,
    "\n")
