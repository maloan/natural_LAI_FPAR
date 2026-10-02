# ==============================================================================
# 00_area_validdomain_after_nonvegetated.R — Create the base valid-domain area
# rasters (0.05° and 0.25°)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(tibble)
  library(readr)
  library(here)
})

source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "cli_args.R"))

terraOptions(progress = 0, memfrac = 0.4)

outdir <- here("analysis", "results", "tables", "masks")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_cdf_if_stale <- function(x, path, dependencies) {
  dependencies <- dependencies[file.exists(dependencies)]
  is_current <- file.exists(path) &&
    length(dependencies) > 0 &&
    file.info(path)$mtime >= max(file.info(dependencies)$mtime)

  if (!is_current) {
    writeCDF(x, path, overwrite = TRUE)
  }
  invisible(path)
}

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
stopifnot(compareGeom(m_nonveg, area, stopOnError = TRUE))

area_sum <- function(cond) {
  global(ifel(cond, area, NA), "sum", na.rm = TRUE)[1, 1] |> as.numeric()
}
support_dom <- is.finite(area) & (area > 0)
nonveg_excl <- support_dom & (m_nonveg == 1)
land_dom <- support_dom & !nonveg_excl & is.finite(m_nonveg)

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

# Write summary and refresh rasters only when their inputs are newer.
write_csv(
  round_numeric(tbl_nonveg, 5),
  file.path(outdir, "domain_nonvegetated_0p05.csv")
)
area_valid_005_file <- here("src", "area_0p05_validdomain_km2.nc")
area_valid_025_file <- here("src", "area_0p25_validdomain_km2.nc")
write_cdf_if_stale(
  area_valid_005,
  area_valid_005_file,
  c(area_005, mask_nonveg)
)
write_cdf_if_stale(
  area_valid_025,
  area_valid_025_file,
  c(area_005, mask_nonveg)
)
