# ==============================================================================
# 02_global_trends_summary.R — Global absolute and relative LAI trends
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(tibble)
  library(here)
})

source(here("R", "helpers", "bootstrap_ci.R"))
source(here("R", "helpers", "cli_args.R"))
source(here("R", "helpers", "io.R"))

vars <- "LAI"
metrics <- c("yearmean", "yearmax")
scenario_spec <- create_scenario_spec(
  cci_alphas = c("alpha_0.05", "alpha_0.1", "alpha_0.2"),
  glc_run_tag = "alpha_0.1"
)
outdir <- here("analysis", "results", "tables", "trends")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

area_file <- here("src", "area_0p25_validdomain_km2.nc")
if (!file.exists(area_file)) {
  stop("Missing area raster: ", area_file)
}
area <- rast(area_file)[[1]]
if (nrow(area) != 720L || ncol(area) != 1440L) {
  stop("Area raster must be a 720 x 1440 global 0.25-degree grid")
}
block_id <- make_block_id(area, block_size_deg = 5)

scenario_weights <- setNames(
  lapply(seq_len(nrow(scenario_spec)), function(i) {
    sc <- scenario_spec[i, ]
    values(
      load_scenario_area(sc$source, sc$run_tag, template = area),
      dataframe = FALSE
    )
  }),
  scenario_spec$scenario
)

summarise_trend_kind <- function(is_relative) {
  results <- list()
  for (var in vars) {
    for (met in metrics) {
      for (i in seq_len(nrow(scenario_spec))) {
        sc <- scenario_spec[i, ]
        trend_values <- read_trend(
          trend_path_factory(
            var,
            met,
            sc$source,
            sc$run_tag,
            is_relative = is_relative
          ),
          sc$scenario,
          template = area
        )
        weights <- scenario_weights[[sc$scenario]]
        if (length(trend_values) != length(weights)) {
          stop("Geometry mismatch for ", sc$scenario)
        }

        valid <- is.finite(trend_values) & is.finite(weights) & weights > 0
        ci <- bootstrap_ci_global(
          x = trend_values[valid],
          w = weights[valid],
          block_id = block_id[valid],
          n_boot = 1000L,
          conf = 0.95
        )
        row <- tibble(
          variable = var,
          metric = met,
          scenario = sc$scenario,
          run_tag = sc$run_tag,
          sig_flag = ci$sig,
          area_km2 = sum(weights[valid]),
          n_pixels = sum(valid)
        )
        if (is_relative) {
          row <- mutate(
            row,
            reltrend_pct_per_year = ci$mean,
            reltrend_ci_lower = ci$lower,
            reltrend_ci_upper = ci$upper,
            reltrend_ci_width = ci$width,
            .before = sig_flag
          )
        } else {
          row <- mutate(
            row,
            abstrend_m2m2yr = ci$mean,
            abstrend_ci_lower = ci$lower,
            abstrend_ci_upper = ci$upper,
            abstrend_ci_width = ci$width,
            .before = sig_flag
          )
        }
        results[[length(results) + 1L]] <- row
      }
    }
  }

  bind_rows(results) |>
    mutate(
      variable = factor(.data$variable, levels = vars),
      metric = factor(.data$metric, levels = metrics),
      scenario = factor(
        .data$scenario,
        levels = scenario_order(scenario_spec)
      )
    ) |>
    arrange(.data$variable, .data$metric, .data$scenario)
}

write_trend_tables <- function(tab, is_relative) {
  kind <- if (is_relative) "relative" else "absolute"
  main <- tab |>
    filter(.data$variable == "LAI", .data$metric == "yearmean") |>
    round_numeric(5)

  write_csv(
    main,
    file.path(
      outdir,
      sprintf("table_global_mean_%s_trends_yearmean_LAI_overview.csv", kind)
    )
  )
  write_csv(
    round_numeric(tab, 5),
    file.path(
      outdir,
      sprintf("global_mean_%s_trends_long_overview.csv", kind)
    )
  )
}

for (is_relative in c(FALSE, TRUE)) {
  write_trend_tables(
    summarise_trend_kind(is_relative),
    is_relative
  )
}
