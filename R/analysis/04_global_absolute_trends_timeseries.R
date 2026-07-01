# ==============================================================================
# 04_global_absolute_trends_timeseries.R — Global annual absolute LAI time
# series and linear trends
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(here)
  library(purrr)
  library(tibble)
})

source(here("R", "helpers", "weighted_means.R"))
source(here("R", "helpers", "plotting.R"))
source(here("R", "helpers", "io.R"))

# Configuration
var <- "LAI"
metrics <- c("yearmean", "yearmax")
year0 <- 1982L
year_end <- 2024L
n_years <- year_end - year0 + 1L
taus <- c("tau_0.05", "tau_0.1", "tau_0.2")
glc_tau <- "tau_0.05"
outdir_fig <- here("analysis", "results", "figures", "timeseries")
outdir_tbl <- here("analysis", "results", "tables", "timeseries")
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tbl, recursive = TRUE, showWarnings = FALSE)

# Area weights
area_path <- here("src", "area_0p25_validdomain_km2.nc")
area <- rast(area_path)[[1]]

make_series <- function(r, area, year0) {
  out <- global_wmean_series(r, area, year0)
  as_tibble(out)
}

compute_ols_full <- function(d,
                             y_col = "value",
                             x_col = "year",
                             robust = TRUE,
                             nw_lag = NULL) {
  d <- d |>
    dplyr::filter(is.finite(.data[[y_col]]), is.finite(.data[[x_col]]))

  if (nrow(d) < 3) {
    return(
      tibble::tibble(
        n = nrow(d),
        slope = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        p_value = NA_real_,
        r2 = NA_real_,
        sig = NA_character_,
        se_type = ifelse(robust, "Newey-West", "OLS"),
        nw_lag = NA_integer_
      )
    )
  }

  m <- stats::lm(stats::formula(paste(y_col, "~", x_col)), data = d)

  compute_trend_ci(
    lm_fit = m,
    robust = robust,
    nw_lag = nw_lag
  )
}

compute_trend_ci <- function(lm_fit,
                             robust = TRUE,
                             nw_lag = NULL) {
  sm <- summary(lm_fit)

  n <- nrow(lm_fit$model)
  slope <- stats::coef(lm_fit)[2]
  r2 <- sm$r.squared

  if (robust) {
    if (is.null(nw_lag)) {
      nw_lag <- floor(4 * (n / 100)^(2 / 9))
    }

    vc <- sandwich::NeweyWest(lm_fit,
      lag = nw_lag,
      prewhite = FALSE,
      adjust = TRUE
    )

    ct <- lmtest::coeftest(lm_fit, vcov. = vc)

    se <- ct[2, 2]
    p_value <- ct[2, 4]
    se_type <- "Newey-West"
  } else {
    se <- sm$coefficients[2, 2]
    p_value <- sm$coefficients[2, 4]
    se_type <- "OLS"
  }

  tcrit <- stats::qt(0.975, df = lm_fit$df.residual)

  tibble::tibble(
    n = n,
    slope = slope,
    ci_lower = slope - tcrit * se,
    ci_upper = slope + tcrit * se,
    p_value = p_value,
    r2 = r2,
    sig = ifelse(p_value < 0.05, "*", ""),
    se_type = se_type,
    nw_lag = ifelse(robust, nw_lag, NA_integer_)
  )
}

predict_trend_line <- function(d, y_col = "value", x_col = "year") {
  d <- d |>
    dplyr::filter(is.finite(.data[[y_col]]), is.finite(.data[[x_col]]))

  if (nrow(d) < 3) {
    return(tibble::tibble(!!x_col := d[[x_col]], fit = NA_real_))
  }

  m <- stats::lm(stats::formula(paste(y_col, "~", x_col)), data = d)

  tibble::tibble(!!x_col := d[[x_col]], fit = stats::predict(m, newdata = d))
}

# Load data
rows <- list()
for (metric in metrics) {
  # Unmasked
  r_unmasked <- load_checked_raster(analysis_raster_path(var, metric, "unmasked", kind = "metric"),
    area,
    n_layers = n_years
  )
  rows[[length(rows) + 1]] <- make_series(r_unmasked, area, year0) |> mutate(metric = metric, scenario = "Unmasked")
  # CCI
  for (tau in taus) {
    r_cci <- load_checked_raster(
      analysis_raster_path(var, metric, "CCI", run_tag = tau, kind = "metric"),
      area,
      n_layers = n_years
    )
    rows[[length(rows) + 1]] <- make_series(r_cci, area, year0) |>
      mutate(metric = metric, scenario = paste0("CCI tau=", sub("^tau_", "", tau)))
  }
  # GLC
  r_glc <- load_checked_raster(
    analysis_raster_path(var, metric, "GLC", run_tag = glc_tau, kind = "metric"),
    area,
    n_layers = n_years
  )
  rows[[length(rows) + 1]] <- make_series(r_glc, area, year0) |> mutate(metric = metric, scenario = "GLC")
}

df <- bind_rows(rows)
df <- df |> mutate(
  metric = factor(
    metric,
    levels = c("yearmean", "yearmax"),
    labels = c("Annual mean", "Annual maximum")
  ),
  scenario = factor(scenario, levels = c(
    "Unmasked", paste0("CCI tau=", sub("^tau_", "", taus)), "GLC"
  ))
)
# OLS trend statistics
trend_stats <- df |>
  group_by(metric, scenario) |>
  group_modify(~ compute_ols_full(.x)) |>
  ungroup()
# Prediction lines
trend_df <- df |>
  nest(data = -c(metric, scenario)) |>
  mutate(trend_pred = map(data, predict_trend_line)) |>
  select(-data) |>
  unnest(trend_pred)
plot_range_df <- df |>
  group_by(metric) |>
  summarise(
    y_max = max(value, na.rm = TRUE),
    y_min = min(value, na.rm = TRUE),
    .groups = "drop"
  )

# Trend annotations
annotation_df <- trend_stats |>
  left_join(plot_range_df, by = "metric") |>
  filter(scenario %in% c("Unmasked", "CCI tau=0.1", "GLC")) |>
  group_by(metric) |>
  arrange(match(scenario, c("Unmasked", "CCI tau=0.1", "GLC")), .by_group = TRUE) |>
  summarise(
    x = 2024.1,
    y = first(y_min) + 0.03 * (first(y_max) - first(y_min)),
    label = {
      lines <- paste0(
        scenario,
        ": slope = ",
        sprintf("%.4f", slope),
        " ± ",
        sprintf("%.4f", (ci_upper - ci_lower) / (2 * stats::qt(0.975, df = n - 2)))
      )
      lines[length(lines)] <- paste0(
        lines[length(lines)],
        "\n(p = ",
        format.pval(p_value[length(p_value)], digits = 4, eps = 1e-5),
        ")"
      )
      paste(lines, collapse = "\n")
    },
    .groups = "drop"
  ) |>
  ungroup()
p <- plot_timeseries(df, trend_df, annotation_df, theme_pub)
# Output
write_csv(
  round_numeric(df, 5),
  file.path(outdir_tbl, "global_timeseries_absolute_trends.csv")
)
write_csv(
  round_numeric(trend_df, 5),
  file.path(outdir_tbl, "global_timeseries_absolute_trends_fit.csv")
)
write_csv(
  round_numeric(trend_stats, 5),
  file.path(
    outdir_tbl,
    "global_timeseries_absolute_trends_statistics.csv"
  )
)
ggsave(
  filename = file.path(outdir_fig, "global_timeseries_absolute_trends.png"),
  plot = p,
  width = 11,
  height = 5.8,
  dpi = 350
)
ggsave(
  filename = file.path(outdir_fig, "global_timeseries_absolute_trends.pdf"),
  plot = p,
  width = 11,
  height = 5.8
)
