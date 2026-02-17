# ==============================================================================
# 01_domain_masking_footprint.R
# Land-use footprints within the valid land domain
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(tibble)
  library(here)
})

# Read inputs
TAU <- "tau_0.2"
OUTDIR <- here("analysis", "results", "masks", TAU)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

tau_num <- as.numeric(sub("^tau_", "", TAU))
tau_token <- sprintf("tau0p%02d", round(tau_num * 100))

MASK_CCI <- here(
  "output", TAU, "masks", "mask_cci",
  sprintf("mask_used_frac_fused_%s_k3_1992-2020_0p05.tif", tau_token)
)

MASK_GLC <- here(
  "output", TAU, "masks", "mask_glc",
  "mask_used_ge3_1992-2020_0p05.tif"
)

MASK_LUH <- here(
  "output", TAU, "masks", "mask_luh_overlap",
  "mask_luh_overlap_CCI_Gmin0p10_Pmin0p10_alpha0p50_1992-2020_0p05_rep.tif"
)


# --- inputs (0.05°) -----------------------------------------------------------
AREA_VALID_005 <- here("src", "area_0p05_validdomain_km2.nc")
A005 <- rast(AREA_VALID_005)[[1]]

m_cci <- rast(MASK_CCI)[[1]] # 0/1 at 0.05
m_glc <- rast(MASK_GLC)[[1]] # 0/1 at 0.05
m_luh <- rast(MASK_LUH)[[1]] # 0/1 at 0.05

compareGeom(m_cci, A005, stopOnError = TRUE)
compareGeom(m_glc, A005, stopOnError = TRUE)
compareGeom(m_luh, A005, stopOnError = TRUE)

# land domain at 0.05, then aggregated domain area per 0.25 cell
land005 <- is.finite(A005) & (A005 > 0)
A_dom_total <- global(ifel(land005, A005, NA), "sum", na.rm = TRUE)[1, 1] |>
  as.numeric()

# --- compute component areas at 0.25° from 0.05° conditions -------------------
summarise_025 <- function(m_lu, label) {
  s <- component_area_global(m_lu, m_luh)
  compact <- tibble::tibble(
    scenario = label,
    denom_land_after_abi_km2 = s$A_dom_total,
    excluded_LU_km2 = s$A_LU,
    excluded_LU_pct = 100 * s$A_LU / s$A_dom_total,
    excluded_LUH_km2 = s$A_LUH,
    excluded_LUH_pct = 100 * s$A_LUH / s$A_dom_total,
    excluded_union_km2 = s$A_union,
    excluded_union_pct = 100 * s$A_union / s$A_dom_total
  )

  partition <- tibble::tibble(
    scenario = label,
    component = c(paste0(
      label, "_only"
    ), "LUH_only", "overlap_LU_and_LUH", "union_any"),
    km2 = c(s$A_LU_only, s$A_LUH_only, s$A_overlap, s$A_union),
    pct_of_land_after_abi = 100 * km2 / s$A_dom_total
  )
  list(compact = compact, partition = partition)
}

component_area_global <- function(LU_005, LUH_005) {
  # make masks safe (binary NA -> 0)
  LU_005[is.na(LU_005)] <- 0
  LUH_005[is.na(LUH_005)] <- 0

  LU <- (LU_005 == 1)
  LUH <- (LUH_005 == 1)

  tot <- function(cond) {
    global(ifel(land005 & cond, A005, NA), "sum", na.rm = TRUE)[1, 1] |>
      as.numeric()
  }

  A_LU_only <- tot(LU & !LUH)
  A_LUH_only <- tot(!LU & LUH)
  A_overlap <- tot(LU & LUH)

  A_LU <- A_LU_only + A_overlap
  A_LUH <- A_LUH_only + A_overlap
  A_union <- A_LU_only + A_LUH_only + A_overlap

  list(
    A_dom_total = A_dom_total,
    A_LU_only = A_LU_only,
    A_LUH_only = A_LUH_only,
    A_overlap = A_overlap,
    A_union = A_union,
    A_LU = A_LU,
    A_LUH = A_LUH
  )
}

cci <- summarise_025(m_cci, "CCI")
glc <- summarise_025(m_glc, "GLC")

tbl_land <- dplyr::bind_rows(cci$compact, glc$compact)
tbl_overlap <- dplyr::bind_rows(cci$partition, glc$partition)

write_csv(
  mutate(tbl_land, across(where(is.double), ~ round(.x, 3))),
  file.path(OUTDIR, sprintf("domain_landuse_after_abiotic_%s_0p25.csv", TAU))
)
write_csv(
  mutate(tbl_overlap, across(where(is.double), ~ round(.x, 3))),
  file.path(OUTDIR, sprintf("domain_overlap_partition_%s_0p25.csv", TAU))
)
