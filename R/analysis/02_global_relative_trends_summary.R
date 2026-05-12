# ==============================================================================
# 02_global_relative_trends_summary.R
# Global mean relative greening
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(tibble)
  library(here)
})

source(here("R", "helpers", "bootstrap_ci.R"))
source(here("R", "helpers", "scenario_config.R"))
source(here("R", "helpers", "io.R"))

# ==============================================================================
# CRS & Coordinate System Validation
# Expected: EPSG:4326 (WGS84 geographic), 0.25° grid (720×1440)
# ==============================================================================

utils::globalVariables(c(
  "variable", "metric", "scenario",
  "reltrend_pct_per_year", "reltrend_ci_lower", "reltrend_ci_upper", "reltrend_ci_width"
))

# Configuration
cci_taus <- c("tau_0.05", "tau_0.1", "tau_0.2")
glc_run_tag <- Sys.getenv("GLC_RUN_TAG", "tau_0.1")

vars <- c("LAI", "FPAR")
metrics <- c("yearmean", "yearmax")

scenario_spec <- create_scenario_spec(cci_taus, glc_run_tag)

dir_unmask <- here("analysis", "unmasked", "0p25")
outdir <- here("analysis", "results", "tables", "trends")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load spatial weights (area in km²)
area_file <- here("src", "area_0p25_validdomain_km2.nc")

if (!file.exists(area_file)) {
  stop("Missing area raster: ", area_file)
}

area <- rast(area_file)[[1]]

# CRS Validation: ensure 0.25° global grid
if (nrow(area) != 720L || ncol(area) != 1440L) {
  stop(
    "Area raster has unexpected dimensions: ", nrow(area), "×", ncol(area),
    "\n  Expected 720×1440 (0.25° grid, EPSG:4326)"
  )
}

area_vals <- values(area, dataframe = FALSE)

area_total <- global(
  ifel(is.finite(area) & area > 0, area, NA),
  "sum",
  na.rm = TRUE
)[1, 1]

if (!is.finite(area_total) || area_total <= 0) {
  stop("Invalid total area.")
}

# Core computation
results <- list()

for (var in vars) {
  for (met in metrics) {
    scenario_results <- list()

    for (i in seq_len(nrow(scenario_spec))) {
      sc <- scenario_spec[i, ]
      path <- trend_path_factory(var, met, sc$source, sc$run_tag, is_relative = TRUE)

      r_vals <- read_trend(path, sc$scenario)

      if (length(r_vals) != length(area_vals)) {
        stop(
          "Geometry mismatch: ", sc$scenario, "\n",
          "  Expected ", length(area_vals), " pixels, got ", length(r_vals)
        )
      }

      ok <- is.finite(r_vals) & is.finite(area_vals) & area_vals > 0

      # Unit conversion: relative trends are fractional (0-1), output as % per year
      r_vals <- r_vals * 100

      # Bootstrap CI for global mean (area-weighted)
      ci <- bootstrap_ci_global(r_vals[ok], area_vals[ok], n_boot = 1000L, conf = 0.95)

      scenario_results[[length(scenario_results) + 1]] <- tibble(
        variable = var,
        metric = met,
        scenario = sc$scenario,
        run_tag = sc$run_tag,
        reltrend_pct_per_year = ci$mean,
        reltrend_ci_lower = ci$lower,
        reltrend_ci_upper = ci$upper,
        reltrend_ci_width = ci$width,
        sig_flag = ci$sig,
        area_km2 = sum(area_vals[ok], na.rm = TRUE),
        n_pixels = sum(ok)
      )
    }

    results[[length(results) + 1]] <- bind_rows(scenario_results)
  }
}

tab <- bind_rows(results) |>
  mutate(
    variable = factor(variable, levels = vars),
    metric = factor(metric, levels = metrics),
    scenario = factor(scenario, levels = scenario_order(scenario_spec))
  ) |>
  arrange(variable, metric, scenario)

# Output formatting and metadata
format_table <- function(df) {
  df |> mutate(across(starts_with("reltrend_"), ~ round(., 4)))
}

# Main table: LAI yearmean
tab_main <- tab |>
  filter(variable == "LAI", metric == "yearmean") |>
  format_table()

write_csv(
  tab_main,
  file.path(outdir, "table_global_mean_relative_trends_yearmean_LAI_overview.csv")
)

# Full overview
write_csv(
  format_table(tab),
  file.path(outdir, "global_mean_relative_trends_long_overview.csv")
)
