suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
cfg <- cfg_read()

run_tag <- cfg$project$run_tag
abi_file <- here(
  "output", run_tag, "masks", "mask_abiotic", "abiotic_components_2007_0p05.tif"
)

r <- rast(abi_file)
area005 <- rast(cfg$grids$grid_005$area_raster)
global_area_km2 <- as.numeric(global(area005, "sum", na.rm = TRUE)[1, 1])

water_drop <- r[[1]] == 1
ice_drop <- r[[2]] == 1

water_only <- water_drop & !ice_drop
ice_only <- ice_drop & !water_drop
overlap <- water_drop & ice_drop
union <- water_drop | ice_drop

area_km2 <- function(mask) {
  as.numeric(global(area005 * mask, "sum", na.rm = TRUE)[1, 1])
}

tbl <- data.frame(
  component = c("water_only", "ice_only", "water_ice_overlap", "union"),
  area_km2 = c(
    area_km2(water_only), area_km2(ice_only), area_km2(overlap), area_km2(union)
  )
)
tbl$area_pct_global <- 100 * tbl$area_km2 / global_area_km2

# rounding
tbl$area_km2 <- round(tbl$area_km2, 0)
tbl$area_pct_global <- round(tbl$area_pct_global, 5)

print(tbl, row.names = FALSE)
