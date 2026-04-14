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

# Read inputs
tau <- Sys.getenv("RUN_TAG", "tau_0.1")
outdir <- here("analysis", "results", "masks", tau)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tau_num <- as.numeric(sub("^tau_", "", tau))
if (!is.finite(tau_num)) {
  stop("RUN_TAG must look like tau_0.1; got: ", tau)
}
tau_token <- sprintf("tau0p%02d", round(tau_num * 100))

pick_latest <- function(dir_path, pattern, label) {
  ff <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (!length(ff)) {
    stop("No ", label, " mask found in ", dir_path, " matching: ", pattern)
  }
  ff[order(file.info(ff)$mtime, decreasing = TRUE)][1]
}

mask_cci <- pick_latest(
  here("output", tau, "masks", "mask_cci"),
  sprintf("^mask_used_frac_fused_%s_.*_0p05\\.tif$", tau_token),
  "CCI"
)

mask_glc <- pick_latest(
  here("output", tau, "masks", "mask_glc"),
  "^mask_used_ge[0-9]+_.*_0p05\\.tif$",
  "GLC"
)

mask_luh <- pick_latest(
  here("output", tau, "masks", "mask_luh_overlap"),
  "^mask_luh_overlap_.*_0p05_rep\\.tif$",
  "LUH overlap"
)

# --- inputs (0.05°) -----------------------------------------------------------
area_valid_005 <- here("src", "area_0p05_validdomain_km2.nc")
if (!file.exists(area_valid_005)) {
  stop("Missing valid-domain area file: ", area_valid_005)
}
a005 <- rast(area_valid_005)[[1]]

m_cci <- rast(mask_cci)[[1]] # 0/1 at 0.05
m_glc <- rast(mask_glc)[[1]] # 0/1 at 0.05
m_luh <- rast(mask_luh)[[1]] # 0/1 at 0.05

compareGeom(m_cci, a005, stopOnError = TRUE)
compareGeom(m_glc, a005, stopOnError = TRUE)
compareGeom(m_luh, a005, stopOnError = TRUE)

# land domain at 0.05, then aggregated domain area per 0.25 cell
land005 <- is.finite(a005) & (a005 > 0)
area_dom_total <- global(ifel(land005, a005, NA), "sum", na.rm = TRUE)[1, 1] |>
  as.numeric()
if (!is.finite(area_dom_total) || area_dom_total <= 0) {
  stop("Invalid domain denominator area (km2): ", area_dom_total)
}

# --- compute component areas at 0.25° from 0.05° conditions -------------------
summarise_025 <- function(m_lu, label) {
  s <- component_area_global(m_lu, m_luh)
  compact <- tibble::tibble(
    scenario = label,
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

cci <- summarise_025(m_cci, "CCI")
glc <- summarise_025(m_glc, "GLC")

tbl_land <- dplyr::bind_rows(cci$compact, glc$compact)
tbl_overlap <- dplyr::bind_rows(cci$partition, glc$partition)

write_csv(
  mutate(tbl_land, across(where(is.double), ~ round(.x, 3))),
  file.path(outdir, sprintf("domain_landuse_after_nonvegetated_%s_0p25.csv", tau))
)
write_csv(
  mutate(tbl_overlap, across(where(is.double), ~ round(.x, 3))),
  file.path(outdir, sprintf("domain_overlap_partition_%s_0p25.csv", tau))
)
