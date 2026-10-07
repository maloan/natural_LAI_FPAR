# =============================================================================
# 06_nonveg_static_from_cci_0p05.R — Build baseline and sensitivity
# non-vegetated masks (0.05°)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})


source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.6)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
cci_dir <- cfg$paths$cci_dir

alpha_water <- cfg$esa_cci$nonvegetated$water_threshold
alpha_ice <- cfg$esa_cci$nonvegetated$ice_threshold
years <- sort(unique(c(1995L, as.integer(cfg$esa_cci$nonvegetated$year), 2022L)))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_nonvegetated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

alpha_w_tok <- gsub("\\.", "p", sprintf("%.2f", alpha_water))
alpha_i_tok <- gsub("\\.", "p", sprintf("%.2f", alpha_ice))

esa_cci <- cfg$esa_cci$classes
vals_water <- as.integer(unlist(esa_cci$water))
vals_ice <- as.integer(unlist(esa_cci$snow_ice))
nodata_vals <- unique(c(as.integer(unlist(esa_cci$nodata)), 255L))

for (year in years) {
  out_tif <- file.path(
    out_dir,
    sprintf(
      "mask_nonvegetated_CCI_%d_alphaW%s_alphaI%s_0p05.tif",
      year,
      alpha_w_tok,
      alpha_i_tok
    )
  )

  in_file <- if (year <= 2015) {
    file.path(
      cci_dir,
      sprintf(
        "ESACCI-LC-L4-LCCS-Map-300m-P1Y-%d-v2.0.7cds.nc",
        year
      )
    )
  } else {
    file.path(
      cci_dir,
      sprintf(
        "C3S-LC-L4-LCCS-Map-300m-P1Y-%d-v2.1.1.nc",
        year
      )
    )
  }

  if (!file.exists(in_file)) {
    stop("Missing expected CCI NetCDF file:\n", in_file)
  }

  components_file <- file.path(out_dir, sprintf("nonvegetated_components_%d_0p05.tif", year))
  dependencies <- c(in_file, here("R", "06_nonveg_static_from_cci_0p05.R"))
  if (file.exists(out_tif) && file.exists(components_file) &&
      min(file.mtime(c(out_tif, components_file))) >= max(file.mtime(dependencies))) {
    message("✓ Non-vegetated mask is current — skipping: ", out_tif)
    next
  }

  message("→ Processing ", basename(in_file), " (year=", year, ")")

  r <- rast(in_file, subds = "lccs_class")
  if (is.na(crs(r))) {
    crs(r) <- crs(ref005)
  }

  r <- terra::subst(r, nodata_vals, NA)

  # Water and ice fractions at 0.05° for this snapshot year.
  pW <- resample(classify(r, cbind(vals_water, 1), others = 0), ref005, method = "average")
  pI <- resample(classify(r, cbind(vals_ice, 1), others = 0), ref005, method = "average")

  water_drop <- ifel(pW >= alpha_water, 1L, 0L)
  ice_drop <- ifel(pI >= alpha_ice, 1L, 0L)
  both_drop <- ifel(water_drop & ice_drop, 1L, 0L)
  nonveg_mask_combined <- ifel(water_drop | ice_drop, 1L, 0L)

  writeRaster(
    c(water_drop, ice_drop, both_drop, nonveg_mask_combined),
    components_file,
    overwrite = TRUE,
    wopt = wopt_byte(FALSE, na = 255L)
  )

  names(nonveg_mask_combined) <- "nonvegetated_drop"
  writeRaster(
    nonveg_mask_combined,
    out_tif,
    overwrite = TRUE,
    wopt = wopt_byte(FALSE, na = 255L)
  )

  rm(r, pW, pI, water_drop, ice_drop, both_drop, nonveg_mask_combined)
  gc(verbose = FALSE)
}
