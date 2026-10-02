# ==============================================================================
# 03_matched_period_annual_mean_trends.R — Unmasked annual-mean LAI trends
# matched to literature comparison periods
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(readr)
  library(tibble)
  library(here)
})

source(here("R", "helpers", "bootstrap_ci.R"))

annual_file <- here(
  "analysis", "unmasked", "0p25", "LAI_georef_yearmean_0p25.nc"
)
area_file <- here("src", "area_0p25_validdomain_km2.nc")
out_file <- here(
  "analysis", "results", "tables", "trends",
  "matched_period_unmasked_annual_mean_LAI_trends.csv"
)

if (!file.exists(annual_file)) {
  stop("Missing annual-mean LAI series: ", annual_file)
}
if (!file.exists(area_file)) {
  stop("Missing post-nonvegetated area raster: ", area_file)
}

# The windows include the periods used for literature comparisons and the
# full-record reference period. All estimates use the same unmasked annual
# mean input and fixed post-nonvegetated areas of valid 0.25-degree cells.
periods <- tribble(
  ~start_year, ~end_year,
  1982L, 2011L,
  1982L, 2015L,
  1982L, 2020L,
  1982L, 2021L,
  2001L, 2017L,
  2001L, 2020L,
  2004L, 2020L,
  1982L, 2024L
)

area <- rast(area_file)[[1]]
area_values <- values(area, dataframe = FALSE)
block_id <- make_block_id(area, block_size_deg = 5)

summarise_period <- function(start_year, end_year) {
  intercept_file <- tempfile(fileext = ".nc")
  slope_file <- tempfile(fileext = ".nc")
  on.exit(unlink(c(intercept_file, slope_file)), add = TRUE)

  status <- system2(
    "cdo",
    c(
      "-O", "trend", sprintf("-selyear,%d/%d", start_year, end_year),
      annual_file, intercept_file, slope_file
    )
  )
  if (status != 0L) {
    stop("CDO trend failed for ", start_year, "--", end_year)
  }

  slope <- rast(slope_file)[[1]]
  compareGeom(slope, area, stopOnError = TRUE)
  trend_values <- values(slope, dataframe = FALSE)
  valid <- is.finite(trend_values) & is.finite(area_values) & area_values > 0
  ci <- bootstrap_ci_global(
    x = trend_values[valid],
    w = area_values[valid],
    block_id = block_id[valid],
    n_boot = 1000L,
    conf = 0.95
  )

  tibble(
    diagnostic = "annual mean",
    start_year = start_year,
    end_year = end_year,
    n_years = end_year - start_year + 1L,
    trend_m2m2yr = ci$mean,
    ci_lower = ci$lower,
    ci_upper = ci$upper,
    n_pixels = sum(valid),
    area_km2 = sum(area_values[valid]),
    n_blocks = ci$n_eff
  )
}

results <- do.call(
  rbind,
  lapply(seq_len(nrow(periods)), function(i) {
    summarise_period(periods$start_year[i], periods$end_year[i])
  })
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write_csv(results, out_file)
message("Wrote ", out_file)
