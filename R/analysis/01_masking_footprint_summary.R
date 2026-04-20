# ==============================================================================
# 01_masking_footprint_summary.R
# Land-use footprints within the valid land domain
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(tibble)
  library(here)
})

source(here("R", "helpers", "files.R"))

# Read inputs
cci_taus <- c("tau_0.05", "tau_0.1", "tau_0.2")
glc_run_tag <- Sys.getenv("GLC_RUN_TAG", "tau_0.1")

outdir <- here("analysis", "results", "masks", "overview")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tau_to_token <- function(run_tag) {
  tau_num <- as.numeric(sub("^tau_", "", run_tag))
  if (!is.finite(tau_num)) {
    stop("run_tag must look like tau_0.1; got: ", run_tag)
  }
  sprintf("tau0p%02d", round(tau_num * 100))
}

load_scenario_masks <- function(run_tag, lu_kind = c("CCI", "GLC")) {
  lu_kind <- match.arg(lu_kind)
  tau_token <- tau_to_token(run_tag)

  lu_dir <- if (lu_kind == "CCI") "mask_cci" else "mask_glc"
  lu_pattern <- if (lu_kind == "CCI") {
    sprintf("^mask_used_frac_fused_%s_.*_0p05\\.tif$", tau_token)
  } else {
    "^mask_used_ge[0-9]+_.*_0p05\\.tif$"
  }

  lu_path <- find_one(
    here("output", run_tag, "masks", lu_dir),
    lu_pattern,
    label = sprintf("%s mask (%s)", lu_kind, run_tag)
  )
  luh_path <- find_one(
    here("output", run_tag, "masks", "mask_luh_overlap"),
    sprintf("^mask_luh_overlap_%s_Gmin0p10_Pmin0p10_alpha0p50_.*_0p05_rep\\.tif$", lu_kind),
    label = sprintf("LUH overlap (%s, %s)", lu_kind, run_tag)
  )

  m_lu <- rast(lu_path)[[1]]
  m_luh <- rast(luh_path)[[1]]
  compareGeom(m_lu, a005, stopOnError = TRUE)
  compareGeom(m_luh, a005, stopOnError = TRUE)

  cat("Selected files for ", lu_kind, " [", run_tag, "]:\n", sep = "")
  cat("  LU:  ", lu_path, "\n", sep = "")
  cat("  LUH: ", luh_path, "\n", sep = "")

  list(lu = m_lu, luh = m_luh)
}

# --- inputs (0.05°) -----------------------------------------------------------
area_valid_005 <- here("src", "area_0p05_validdomain_km2.nc")
if (!file.exists(area_valid_005)) {
  stop("Missing valid-domain area file: ", area_valid_005)
}
a005 <- rast(area_valid_005)[[1]]
land005 <- is.finite(a005) & (a005 > 0)
area_dom_total <- global(ifel(land005, a005, NA), "sum", na.rm = TRUE)[1, 1] |>
  as.numeric()
if (!is.finite(area_dom_total) || area_dom_total <= 0) {
  stop("Invalid domain denominator area (km2): ", area_dom_total)
}

# --- compute component areas at 0.25° from 0.05° conditions -------------------
summarise_025 <- function(m_lu, m_luh, label, run_tag) {
  s <- component_area_global(m_lu, m_luh)
  compact <- tibble::tibble(
    scenario = label,
    run_tag = run_tag,
    denom_land_after_abi_km2 = s$area_dom_total,
    excluded_LU_km2 = s$area_lu,
    excluded_LU_pct = 100 * s$area_lu / s$area_dom_total,
    excluded_LUH_km2 = s$area_luh,
    excluded_LUH_pct = 100 * s$area_luh / s$area_dom_total,
    excluded_union_km2 = s$area_union,
    excluded_union_pct = 100 * s$area_union / s$area_dom_total
  )

  partition <- tibble::tibble(
    scenario = label,
    run_tag = run_tag,
    component = c(paste0(
      label, "_only"
    ), "LUH_only", "overlap_LU_and_LUH", "union_any"),
    km2 = c(s$area_lu_only, s$area_luh_only, s$area_overlap, s$area_union),
    pct_of_land_after_abi = 100 * .data$km2 / s$area_dom_total
  )
  list(compact = compact, partition = partition)
}

component_area_global <- function(lu_005, luh_005) {
  # make masks safe (binary NA -> 0)
  lu_005[is.na(lu_005)] <- 0
  luh_005[is.na(luh_005)] <- 0

  lu <- (lu_005 == 1)
  luh <- (luh_005 == 1)

  tot <- function(cond) {
    global(ifel(land005 & cond, a005, NA), "sum", na.rm = TRUE)[1, 1] |>
      as.numeric()
  }

  area_lu_only <- tot(lu & !luh)
  area_luh_only <- tot(!lu & luh)
  area_overlap <- tot(lu & luh)

  area_lu <- area_lu_only + area_overlap
  area_luh <- area_luh_only + area_overlap
  area_union <- area_lu_only + area_luh_only + area_overlap

  list(
    area_dom_total = area_dom_total,
    area_lu_only = area_lu_only,
    area_luh_only = area_luh_only,
    area_overlap = area_overlap,
    area_union = area_union,
    area_lu = area_lu,
    area_luh = area_luh
  )
}

all_compact <- list()
all_partition <- list()

for (tau in cci_taus) {
  ms <- load_scenario_masks(run_tag = tau, lu_kind = "CCI")
  lab <- sprintf("CCI tau=%s", sub("^tau_", "", tau))
  s <- summarise_025(ms$lu, ms$luh, lab, tau)
  all_compact[[length(all_compact) + 1]] <- s$compact
  all_partition[[length(all_partition) + 1]] <- s$partition
}

ms_glc <- load_scenario_masks(run_tag = glc_run_tag, lu_kind = "GLC")
s_glc <- summarise_025(ms_glc$lu, ms_glc$luh, "GLC", glc_run_tag)
all_compact[[length(all_compact) + 1]] <- s_glc$compact
all_partition[[length(all_partition) + 1]] <- s_glc$partition

tbl_land <- dplyr::bind_rows(all_compact)
tbl_overlap <- dplyr::bind_rows(all_partition)

scenario_order <- c(
  "CCI tau=0.05",
  "CCI tau=0.1",
  "CCI tau=0.2",
  "GLC"
)

tbl_land <- tbl_land |>
  mutate(
    scenario = factor(scenario, levels = scenario_order)
  ) |>
  arrange(scenario, run_tag) |>
  mutate(scenario = as.character(scenario))

tbl_overlap <- tbl_overlap |>
  mutate(
    scenario = factor(scenario, levels = scenario_order)
  ) |>
  arrange(scenario, run_tag, component) |>
  mutate(scenario = as.character(scenario))

write_csv(
  mutate(tbl_land, across(where(is.double), ~ round(.x, 3))),
  file.path(outdir, "domain_landuse_after_nonvegetated_overview_0p25.csv")
)
write_csv(
  mutate(tbl_overlap, across(where(is.double), ~ round(.x, 3))),
  file.path(outdir, "domain_overlap_partition_overview_0p25.csv")
)
