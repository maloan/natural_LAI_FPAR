# ==============================================================================
# 02_global_mean_relative_trends.R
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
TAU <- "tau_0.1"

DIR_UNMASK <- here("analysis", "unmasked", "0p25")
DIR_EVAL <- here("output", TAU, "eval")

OUTDIR <- here("analysis", "results", "global_mean_relative_trends", TAU)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

VARS <- c("LAI", "FPAR")
METRICS <- c("yearmean", "yearmax")

case_labels <- c(
  unmasked   = "Unmasked (post-abiotic)",
  masked_CCI = "Natural-only (CCI-based)",
  masked_GLC = "Natural-only (GLC-based)"
)

# ---- valid-domain area raster --------------------------------
AREA_025 <- here("src", "area_0p25_validdomain_km2.nc")
A <- rast(AREA_025)[[1]]

# ---- helpers -----------------------------------------------------------------
wmean_global_ok <- function(x, ok) {
  if (is.null(x)) {
    return(NA_real_)
  }
  num <- global(ifel(ok, x * A, NA), "sum", na.rm = TRUE)[1, 1]
  den <- global(ifel(ok, A, NA), "sum", na.rm = TRUE)[1, 1]
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  as.numeric(num / den)
}

area_ok_km2 <- function(ok) {
  global(ifel(ok, A, NA), "sum", na.rm = TRUE)[1, 1] |>
    as.numeric()
}
round_df <- function(df) mutate(df, across(where(is.double), ~ round(.x, 3)))

mk_table_yearmean <- function(tab, var_keep) {
  tab |>
    filter(variable == var_keep, metric == "yearmean") |>
    transmute(
      Variable = as.character(variable),
      Metric = as.character(metric),
      Domain = as.character(case_label),
      `Effective area (million km²)` = area_km2 / 1e6,
      `Global mean relative trend (% yr^-1)` = reltrend_pct_per_year,
      effective_area_pct = 100 * area_km2 / 119284029
    )
}

# ---- compute long table ------------------------------------------------------
rows <- list()

for (var in VARS) {
  for (met in METRICS) {
    rU <- rast(file.path(
      DIR_UNMASK,
      sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, met)
    ))
    rC <- rast(file.path(
      DIR_EVAL, sprintf("trend_%s_%s", var, "CCI"),
      sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, met)
    ))
    rG <- rast(file.path(
      DIR_EVAL, sprintf("trend_%s_%s", var, "GLC"),
      sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, met)
    ))
    compareGeom(A, rU, stopOnError = TRUE)
    compareGeom(A, rC, stopOnError = TRUE)
    compareGeom(A, rG, stopOnError = TRUE)

    okU <- (is.finite(A) & A > 0) & is.finite(rU)
    okC <- (is.finite(A) & A > 0) & is.finite(rC)
    okG <- (is.finite(A) & A > 0) & is.finite(rG)

    rows[[length(rows) + 1]] <- tibble(
      variable = var,
      metric = met,
      case = c("unmasked", "masked_CCI", "masked_GLC"),
      reltrend_pct_per_year = 100 * c(
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
    variable = factor(variable, levels = VARS),
    metric = factor(metric, levels = METRICS),
    case = factor(case, levels = names(case_labels)),
    case_label = factor(
      case_labels[as.character(case)],
      levels = unname(case_labels)
    )
  ) |>
  arrange(variable, metric, case_label)

# ---- tables  --------------------------------
tab_main <- round_df(mk_table_yearmean(tab, "LAI"))

write_csv(
  tab_main, file.path(
    OUTDIR, sprintf(
      "table_global_mean_relative_trends_yearmean_LAI_%s.csv", TAU
    )
  )
)
write_csv(
  round_df(tab), file.path(
    OUTDIR, sprintf("global_mean_relative_trends_long_%s.csv", TAU)
  )
)
