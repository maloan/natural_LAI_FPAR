# ==============================================================================
# 03_global_absolute_trends_summary.R
# Global mean absolute greening
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(tibble)
  library(here)
})

# ---- config ------------------------------------------------------------------
tau <- Sys.getenv("RUN_TAG", "tau_0.05")

dir_unmask <- here("analysis", "unmasked", "0p25")
dir_eval <- here("output", tau, "eval")

outdir <- here("analysis", "results", "global_mean_absolute_trends", tau)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

vars <- c("LAI", "FPAR")
metrics <- c("yearmean", "yearmax")

case_labels <- c(
  unmasked   = "Unmasked (post-non-vegetated)",
  masked_CCI = "Natural-only (CCI-based)",
  masked_GLC = "Natural-only (GLC-based)"
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

read_single_layer_trend <- function(path, label) {
  r <- rast(path)
  if (nlyr(r) != 1) {
    stop(
      sprintf("Expected exactly 1 layer in %s, got %d: %s", label, nlyr(r), path)
    )
  }
  r[[1]]
}

mk_table_yearmean <- function(tab, var_keep) {
  out <- tab |>
    filter(.data$variable == var_keep, .data$metric == "yearmean") |>
    transmute(
      Variable = as.character(.data$variable),
      Metric = as.character(.data$metric),
      Domain = as.character(.data$case_label),
      `Effective area (million km²)` = .data$area_km2 / 1e6,
      trend_abs_per_year = .data$abs_trend_per_year,
        effective_area_pct = 100 * .data$area_km2 / area_dom_total
    )
  names(out)[names(out) == "trend_abs_per_year"] <- sprintf(
    "Global mean absolute trend (%s yr^-1)",
    var_keep
  )
  out
}

# ---- compute long table ------------------------------------------------------
rows <- list()

for (var in vars) {
  for (met in metrics) {
    f_U <- file.path(
      dir_unmask,
      sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, met)
    )
    f_C <- file.path(
      dir_eval, sprintf("trend_%s_%s", var, "CCI"),
      sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, met)
    )
    f_G <- file.path(
      dir_eval, sprintf("trend_%s_%s", var, "GLC"),
      sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, met)
    )

    if (!file.exists(f_U)) stop("Unmasked trend file missing: ", f_U)
    if (!file.exists(f_C)) stop("CCI trend file missing: ", f_C)
    if (!file.exists(f_G)) stop("GLC trend file missing: ", f_G)

    rU <- read_single_layer_trend(f_U, "unmasked trend")
    rC <- read_single_layer_trend(f_C, "CCI trend")
    rG <- read_single_layer_trend(f_G, "GLC trend")
    compareGeom(area, rU, stopOnError = TRUE)
    compareGeom(area, rC, stopOnError = TRUE)
    compareGeom(area, rG, stopOnError = TRUE)

    okU <- (is.finite(area) & area > 0) & is.finite(rU)
    okC <- (is.finite(area) & area > 0) & is.finite(rC)
    okG <- (is.finite(area) & area > 0) & is.finite(rG)

    rows[[length(rows) + 1]] <- tibble(
      variable = var,
      metric = met,
      case = c("unmasked", "masked_CCI", "masked_GLC"),
      abs_trend_per_year = c(
        wmean_global_ok(rU, okU),
        wmean_global_ok(rC, okC),
        wmean_global_ok(rG, okG)
      ),
      area_km2 = c(
        area_ok_km2(okU),
        area_ok_km2(okC),
        area_ok_km2(okG)
      )
    )
  }
}

tab <- bind_rows(rows) |>
  mutate(
    variable = factor(variable, levels = vars),
    metric = factor(metric, levels = metrics),
    case = factor(case, levels = names(case_labels)),
    case_label = factor(
      case_labels[as.character(case)],
      levels = unname(case_labels)
    )
  ) |>
  arrange(variable, metric, case_label)

# ---- tables  --------------------------------
for (var in vars) {
  tab_main <- round_df(mk_table_yearmean(tab, var))
  write_csv(
    tab_main, file.path(
      outdir, sprintf(
        "table_global_mean_absolute_trends_yearmean_%s_%s.csv", var, tau
      )
    )
  )
}

write_csv(
  round_df(tab), file.path(
    outdir, sprintf("global_mean_absolute_trends_long_%s.csv", tau)
  )
)
