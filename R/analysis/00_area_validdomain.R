# ==============================================================================
# 00_area_validdomain_from_abiotic.R
# Create valid-domain area rasters after abiotic masking (0.05° and 0.25°)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(tibble)
  library(readr)
  library(here)
})

OUTDIR <- here("analysis", "results", "masks")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

AREA_005 <- here("src", "area_0p05_km2.nc")

MASK_ABIOTIC <- here(
  "output", "tau_0.1", "masks", "mask_abiotic",
  "mask_abiotic_CCI_2007_tauW0p05_tauI0p05_0p05.tif"
)

A <- rast(AREA_005)[[1]]
m_abi <- rast(MASK_ABIOTIC)[[1]]
compareGeom(m_abi, A, stopOnError = TRUE)

area_sum <- function(cond) {
  global(ifel(cond, A, NA), "sum", na.rm = TRUE)[1, 1] |> as.numeric()
}

support_dom <- is.finite(A) & (A > 0)
m_abi[is.na(m_abi)] <- 0
abi_excl <- (m_abi == 1)
land_dom <- support_dom & !abi_excl

A_support <- area_sum(support_dom)
A_abi_excl <- area_sum(support_dom & abi_excl)
A_land_after_abi <- area_sum(land_dom)

tbl_abi <- tibble(
  denom_support_km2 = A_support,
  abiotic_removed_km2 = A_abi_excl,
  abiotic_removed_pct_of_support = 100 * (A_abi_excl / A_support),
  land_after_abi_km2 = A_land_after_abi,
  land_after_abi_pct_of_support = 100 * (A_land_after_abi / A_support)
)

write_csv(tbl_abi, file.path(OUTDIR, "domain_abiotic_0p05.csv"))

area_valid_005 <- A
area_valid_005[!land_dom] <- NA

writeCDF(area_valid_005,
  here("src", "area_0p05_validdomain_km2.nc"),
  overwrite = TRUE
)

area_valid_025 <- aggregate(area_valid_005, fact = 5, fun = "sum", na.rm = TRUE)
area_valid_025[area_valid_025 == 0] <- NA

writeCDF(area_valid_025,
  here("src", "area_0p25_validdomain_km2.nc"),
  overwrite = TRUE
)
