# =============================================================================
# 04_glc_stack_0p05.R — Build annual GLC_FCS30D categorical yearstack (0.05°)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "plotting.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.9)

ref005 <- rast(cfg$grids$grid_005$ref_raster)

glc_dir <- cfg$paths$glc_dir
out_dir <- cfg$paths$glc_out_dir
stack_out <- file.path(out_dir, "glc_cat_yearstack_0p05.tif")
grass_out <- file.path(out_dir, "glc_grass_yearstack_0p05.tif")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nodata_vals <- as.integer(unlist(cfg$glc$classes$nodata))

years <- as.integer(cfg$glc$years)
annual_files <- file.path(
  glc_dir,
  sprintf("GLC_FCS30D_mode_0p05_%d.tif", years)
)
missing_files <- annual_files[!file.exists(annual_files)]
if (length(missing_files)) {
  stop(
    "Missing expected GLC annual files:\n",
    paste(missing_files, collapse = "\n")
  )
}

#  rebuild stack
bands <- vector("list", length(years))

for (i in seq_along(years)) {
  yr <- years[i]
  annual_file <- annual_files[i]

  message("→ Processing ", basename(annual_file))
  r <- rast(annual_file)[[1]]

  if (is.na(crs(r))) {
    crs(r) <- crs(ref005)
  }

  r <- align_to_template(r, ref005, method = "near")

  if (length(nodata_vals)) {
    r <- terra::subst(r, nodata_vals, NA)
  }

  names(r) <- sprintf("Y%04d", yr)
  bands[[i]] <- r

  gc()
}

stack <- rast(bands)
writeRaster(stack, stack_out, overwrite = TRUE)

grass_file <- file.path(
  glc_dir,
  "GLC_FCS30D_grass_fraction_0p05_1985_2022.tif"
)
if (!file.exists(grass_file)) {
  stop("Missing expected GLC grass-fraction file:\n", grass_file)
}

message("→ Processing ", basename(grass_file))
grass <- rast(grass_file)

if (nlyr(grass) != length(years)) {
  stop(
    "Expected ", length(years), " grass-fraction layers, found ",
    nlyr(grass)
  )
}

grass <- align_to_template(grass, ref005, method = "bilinear")
grass <- clamp(grass, 0, 1)
names(grass) <- sprintf("Y%04d", years)
writeRaster(grass, grass_out, overwrite = TRUE)
gc()
