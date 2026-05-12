# ==============================================================================
# 04_global_absolute_trends_timeseries.R
# Global annual absolute LAI time series and linear trends
#
# Temporal coverage: 1982–2024 (43 years)
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
source(here("R", "helpers", "statistics.R"))

# ==============================================================================
# Configuration
# ==============================================================================

var <- "LAI"

metrics <- c(
  "yearmean",
  "yearmax"
)

year0 <- 1982L
year_end <- 2024L
n_years <- year_end - year0 + 1L

taus <- c(
  "tau_0.05",
  "tau_0.1",
  "tau_0.2"
)

# Note: GLC uses only tau_0.1 as the primary masking scenario, while CCI
# includes all three tau values to explore sensitivity to the non-vegetated area
# threshold.
glc_tau <- "tau_0.1"

outdir <- here(
  "analysis",
  "results",
  "figures",
  "timeseries"
)

if (!dir.create(outdir, recursive = TRUE, showWarnings = FALSE)) {
  if (!dir.exists(outdir)) {
    stop(
      sprintf("Failed to create output directory: %s", outdir),
      call. = FALSE
    )
  }
}

# ==============================================================================
# Area weights
# ==============================================================================

area_path <- here(
  "src",
  "area_0p25_validdomain_km2.nc"
)

if (!file.exists(area_path)) {
  stop("Missing area raster: ", area_path, call. = FALSE)
}

area <- rast(area_path)[[1]]

# ==============================================================================
# Validation helpers
# ==============================================================================

check_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path, call. = FALSE)
  }

  path
}

load_checked <- function(path, area, n_years) {
  path <- check_file(path)

  r <- rast(path)

  if (nlyr(r) != n_years) {
    stop(
      sprintf(
        "Expected %d layers, found %d in %s",
        n_years,
        nlyr(r),
        basename(path)
      ),
      call. = FALSE
    )
  }

  compareGeom(r[[1]], area, stopOnError = TRUE)

  r
}

make_series <- function(r, area, year0) {
  out <- global_wmean_series(r, area, year0)

  if (nrow(out) == 0) {
    stop("Weighted mean time series is empty", call. = FALSE)
  }

  as_tibble(out)
}

# ==============================================================================
# File paths
# ==============================================================================

path_unmasked <- function(metric) {
  here(
    "analysis",
    "unmasked",
    "0p25",
    sprintf(
      "%s_georef_%s_0p25.nc",
      var,
      metric
    )
  )
}

path_cci <- function(metric, tau) {
  here(
    "output",
    tau,
    "eval",
    sprintf("trend_%s_CCI", var),
    sprintf(
      "%s_%s_0p25.nc",
      var,
      metric
    )
  )
}

path_glc <- function(metric, tau) {
  here(
    "output",
    tau,
    "eval",
    sprintf("trend_%s_GLC", var),
    sprintf(
      "%s_%s_0p25.nc",
      var,
      metric
    )
  )
}

# ==============================================================================
# Load data
# ==============================================================================

rows <- list()

for (metric in metrics) {
  # ---------------------------------------------------------------------------
  # Unmasked
  # ---------------------------------------------------------------------------

  r_unmasked <- load_checked(
    path_unmasked(metric),
    area,
    n_years
  )

  rows[[length(rows) + 1]] <-
    make_series(r_unmasked, area, year0) |>
    mutate(
      metric = metric,
      scenario = "Unmasked"
    )

  # ---------------------------------------------------------------------------
  # CCI
  # ---------------------------------------------------------------------------

  for (tau in taus) {
    r_cci <- load_checked(
      path_cci(metric, tau),
      area,
      n_years
    )

    rows[[length(rows) + 1]] <-
      make_series(r_cci, area, year0) |>
      mutate(
        metric = metric,
        scenario = paste0(
          "CCI tau=",
          sub("^tau_", "", tau)
        )
      )
  }

  # ---------------------------------------------------------------------------
  # GLC
  # ---------------------------------------------------------------------------

  r_glc <- load_checked(
    path_glc(metric, glc_tau),
    area,
    n_years
  )

  rows[[length(rows) + 1]] <-
    make_series(r_glc, area, year0) |>
    mutate(
      metric = metric,
      scenario = "GLC"
    )
}

df <- bind_rows(rows)

if (nrow(df) == 0) {
  stop("No data loaded from files", call. = FALSE)
}

df <- df |>
  mutate(
    metric = factor(
      metric,
      levels = c(
        "yearmean",
        "yearmax"
      ),
      labels = c(
        "Annual mean",
        "Annual maximum"
      )
    ),
    scenario = factor(
      scenario,
      levels = c(
        "Unmasked",
        paste0(
          "CCI tau=",
          sub("^tau_", "", taus)
        ),
        "GLC"
      )
    )
  )

# ==============================================================================
# OLS trend statistics
# ==============================================================================

trend_stats <- df |>
  group_split(metric, scenario) |>
  map_dfr(compute_ols_full)

# ==============================================================================
# Prediction lines
# ==============================================================================

trend_df <- df |>
  nest(data = -c(metric, scenario)) |>
  mutate(trend_pred = map(data, predict_trend_line)) |>
  select(-data) |>
  unnest(trend_pred)

# ==============================================================================
# Plot
# ==============================================================================

# Define color palette
cols <- c(
  "Unmasked" = "grey20",
  "CCI tau=0.05" = "#1b9e77",
  "CCI tau=0.1" = "#d95f02",
  "CCI tau=0.2" = "#7570b3",
  "GLC" = "#386cb0"
)

# Verify all scenario levels have colors defined
scenario_levels <- levels(df$scenario)
missing_colors <- setdiff(scenario_levels, names(cols))
if (length(missing_colors) > 0) {
  stop(
    sprintf(
      "Color palette missing entries for: %s",
      paste(missing_colors, collapse = ", ")
    ),
    call. = FALSE
  )
}

p <- ggplot() +
  geom_ribbon(
    data = trend_df,
    aes(
      x = year,
      ymin = ci_lower,
      ymax = ci_upper,
      fill = scenario
    ),
    alpha = 0.15
  ) +
  geom_line(
    data = df,
    aes(
      x = year,
      y = value,
      colour = scenario
    ),
    alpha = 0.35,
    linewidth = 0.3
  ) +
  geom_line(
    data = trend_df,
    aes(
      x = year,
      y = fit,
      colour = scenario
    ),
    linewidth = 0.9
  ) +
  facet_wrap(
    ~metric,
    scales = "fixed"
  ) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(
    values = cols,
    guide = "none"
  ) +
  labs(
    x = "Year",
    y = "Global mean LAI (dimensionless)",
    title = "Global LAI trajectories and linear trends (1982–2024)",
    subtitle = "Shaded bands: 95% confidence interval of fitted mean response"
  ) +
  theme_pub(base_size = 11)

# ==============================================================================
# Output
# ==============================================================================

write_csv(
  df,
  file.path(
    outdir,
    "global_timeseries_absolute_trends.csv"
  )
)

write_csv(
  trend_df,
  file.path(
    outdir,
    "global_timeseries_absolute_trends_fit.csv"
  )
)

write_csv(
  trend_stats,
  file.path(
    outdir,
    "global_timeseries_absolute_trends_statistics.csv"
  )
)

ggsave(
  filename = file.path(
    outdir,
    "global_timeseries_absolute_trends.png"
  ),
  plot = p,
  width = 11,
  height = 5.8,
  dpi = 350
)

ggsave(
  filename = file.path(
    outdir,
    "global_timeseries_absolute_trends.pdf"
  ),
  plot = p,
  width = 11,
  height = 5.8
)
