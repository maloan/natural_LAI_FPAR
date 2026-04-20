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

# ---- config ------------------------------------------------------------------
cci_taus <- c("tau_0.05", "tau_0.1", "tau_0.2")
glc_run_tag <- Sys.getenv("GLC_RUN_TAG", "tau_0.1")

dir_unmask <- here("analysis", "unmasked", "0p25")

outdir <- here("analysis", "results", "global_mean_relative_trends", "overview")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

vars <- c("LAI", "FPAR")
metrics <- c("yearmean", "yearmax")
scenario_order <- c("Unmasked", "CCI tau=0.05", "CCI tau=0.1", "CCI tau=0.2", "GLC")

scenario_spec <- tibble(
  scenario = scenario_order,
  source = c("unmasked", "CCI", "CCI", "CCI", "GLC"),
  run_tag = c(NA_character_, cci_taus, glc_run_tag)
)

# ---- valid-domain area raster --------------------------------
area_025 <- here("src", "area_0p25_validdomain_km2.nc")
if (!file.exists(area_025)) {
  stop("Missing valid-domain area raster: ", area_025)
}
area <- rast(area_025)[[1]]
area_dom_total <- global(ifel(is.finite(area) & area > 0, area, NA), "sum", na.rm = TRUE)[1, 1] |>
  as.numeric()
if (!is.finite(area_dom_total) || area_dom_total <= 0) {
  stop("Invalid valid-domain area denominator (km2): ", area_dom_total)
}

# ---- helpers -----------------------------------------------------------------
wmean_global_ok <- function(x, ok) {
  if (is.null(x)) {
    return(NA_real_)
  }
  num <- global(ifel(ok, x * area, NA), "sum", na.rm = TRUE)[1, 1]
  den <- global(ifel(ok, area, NA), "sum", na.rm = TRUE)[1, 1]
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  as.numeric(num / den)
}

area_ok_km2 <- function(ok) {
  global(ifel(ok, area, NA), "sum", na.rm = TRUE)[1, 1] |>
    as.numeric()
}
round_df <- function(df) mutate(df, across(where(is.double), ~ round(.x, 3)))

mk_table_yearmean <- function(tab, var_keep) {
  tab |>
    filter(.data$variable == var_keep, .data$metric == "yearmean") |>
    transmute(
      Variable = as.character(.data$variable),
      Metric = as.character(.data$metric),
      Domain = as.character(.data$scenario),
      `Run tag` = as.character(.data$run_tag),
      `Effective area (million km²)` = .data$area_km2 / 1e6,
      `Global mean relative trend (% yr^-1)` = .data$reltrend_pct_per_year,
        effective_area_pct = 100 * .data$area_km2 / area_dom_total
    )
}

trend_path <- function(var, met, source, run_tag = NULL) {
  if (identical(source, "unmasked")) {
    file.path(
      dir_unmask,
      sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, met)
    )
  } else {
    file.path(
      here("output", run_tag, "eval", sprintf("trend_%s_%s", var, source)),
      sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, met)
    )
  }
}

# ---- compute long table ------------------------------------------------------
rows <- list()

for (var in vars) {
  for (met in metrics) {
    scenario_rows <- list()
    for (i in seq_len(nrow(scenario_spec))) {
      sc <- scenario_spec[i, ]
      f_path <- trend_path(var, met, sc$source, sc$run_tag)
      if (!file.exists(f_path)) {
        stop("Missing trend file for ", sc$scenario, ": ", f_path)
      }

      r <- rast(f_path)
      compareGeom(area, r, stopOnError = TRUE)

      ok <- (is.finite(area) & area > 0) & is.finite(r)
      scenario_rows[[length(scenario_rows) + 1]] <- tibble(
        variable = var,
        metric = met,
        scenario = sc$scenario,
        run_tag = sc$run_tag,
        reltrend_pct_per_year = 100 * wmean_global_ok(r, ok),
        area_km2 = area_ok_km2(ok)
      )
    }

    rows[[length(rows) + 1]] <- bind_rows(scenario_rows)
  }
}

tab <- bind_rows(rows) |>
  mutate(
    variable = factor(variable, levels = vars),
    metric = factor(metric, levels = metrics),
    scenario = factor(scenario, levels = scenario_order)
  ) |>
  arrange(variable, metric, scenario)

# ---- tables  --------------------------------
tab_main <- round_df(mk_table_yearmean(tab, "LAI"))

write_csv(
  tab_main, file.path(
    outdir, "table_global_mean_relative_trends_yearmean_LAI_overview.csv"
  )
)
write_csv(
  round_df(tab), file.path(
    outdir, "global_mean_relative_trends_long_overview.csv"
  )
)
