# ==============================================================================
# Area valid-domain summary after non-vegetated masking
# Create valid-domain area rasters after non-vegetated masking (0.05° and 0.25°)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(tibble)
  library(readr)
  library(here)
})

outdir <- here("analysis", "results", "masks")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

area_005 <- here("src", "area_0p05_km2.nc")
run_tag <- Sys.getenv("RUN_TAG", "tau_0.1")

mask_dir <- here("output", run_tag, "masks", "mask_nonvegetated")
mask_candidates <- list.files(
  mask_dir,
  pattern = "mask_nonvegetated_CCI_2007_.*_0p05\\.tif$",
  full.names = TRUE
)
if (!length(mask_candidates)) {
  stop("No nonvegetated mask found in: ", mask_dir)
}
mask_nonveg <- mask_candidates[order(file.info(mask_candidates)$mtime, decreasing = TRUE)][1]

if (!file.exists(area_005)) {
  stop("Missing area raster: ", area_005)
}

area <- rast(area_005)[[1]]
m_nonveg <- rast(mask_nonveg)[[1]]
compareGeom(m_nonveg, area, stopOnError = TRUE)

area_sum <- function(cond) {
  global(ifel(cond, area, NA), "sum", na.rm = TRUE)[1, 1] |> as.numeric()
}

support_dom <- is.finite(area) & (area > 0)
m_nonveg[is.na(m_nonveg)] <- 0
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

write_csv(tbl_nonveg, file.path(outdir, "domain_nonvegetated_0p05.csv"))

area_valid_005 <- area
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
