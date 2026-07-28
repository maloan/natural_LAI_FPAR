# ==============================================================================
# 00_area_validdomain_after_nonvegetated.R — Create valid-domain area rasters
# after non-vegetated masking (0.05° and 0.25°)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(tibble)
  library(readr)
  library(here)
})

source(here("R", "helpers", "io.R"))

outdir <- here("analysis", "results", "tables", "masks")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

area_005 <- here("src", "area_0p05_km2.nc")
mask_nonveg <- here(
  "output",
  "alpha_0.1",
  "masks",
  "mask_nonvegetated",
  "mask_nonvegetated_CCI_2007_alphaW0p05_alphaI0p05_0p05.tif"
)

area <- rast(area_005)[[1]]
m_nonveg <- rast(mask_nonveg)[[1]]
compareGeom(m_nonveg, area, stopOnError = TRUE)

area_sum <- function(cond) {
  global(ifel(cond, area, NA), "sum", na.rm = TRUE)[1, 1] |> as.numeric()
}
support_dom <- is.finite(area) & (area > 0)
nonveg_excl <- support_dom & (m_nonveg == 1)
land_dom <- support_dom & !nonveg_excl & is.finite(m_nonveg)
nonveg_excl <- (m_nonveg == 1)
land_dom <- support_dom & !nonveg_excl

area_support <- area_sum(support_dom)
area_nonveg_excl <- area_sum(support_dom & nonveg_excl)
area_land_after_nonveg <- area_sum(land_dom)

tbl_nonveg <- tibble(
  denom_support_km2 = area_support,
  nonvegetated_removed_km2 = area_nonveg_excl,
  nonvegetated_removed_pct_of_support = 100 * (area_nonveg_excl / area_support),
  land_after_nonvegetated_km2 = area_land_after_nonveg,
  land_after_nonvegetated_pct_of_support = 100 * (area_land_after_nonveg / area_support)
)

area_valid_005 <- area
area_valid_005[!land_dom] <- NA
area_valid_025 <- aggregate(area_valid_005,
  fact = 5,
  fun = "sum",
  na.rm = TRUE
)
area_valid_025[area_valid_025 == 0] <- NA

# Write files
write_csv(
  round_numeric(tbl_nonveg, 5),
  file.path(outdir, "domain_nonvegetated_0p05.csv")
)
writeCDF(area_valid_005,
  here("src", "area_0p05_validdomain_km2.nc"),
  overwrite = TRUE
)
writeCDF(area_valid_025,
  here("src", "area_0p25_validdomain_km2.nc"),
  overwrite = TRUE
)
