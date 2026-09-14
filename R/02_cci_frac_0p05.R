## =============================================================================
# 02_cci_frac_0p05.R — Aggregate ESA-CCI/C3S land cover to 0.05° fractional cover
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.6)

cci_dir <- cfg$paths$cci_dir
out_dir <- cfg$paths$cci_out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tmpl <- rast(cfg$grids$grid_005$ref_raster)

esa_cci <- cfg$esa_cci$classes
nodata_vals <- unique(c(esa_cci$nodata, 255))

groups <- Filter(
  Negate(is.null),
  list(
    cropland  = esa_cci$cropland,
    urban     = esa_cci$urban,
    cls30     = esa_cci$cls30,
    cls40     = esa_cci$cls40,
    grassland = esa_cci$grassland
  )
)

start_year <- cfg$project$years$cci_start
end_year <- cfg$project$years$cci_end

# Precompute write options and GDAL settings once
gdal_opts <- gdal_co_f32(FALSE)

plan_year <- start_year:end_year
plan_name <- ifelse(
  plan_year <= 2015,
  sprintf(
    "ESACCI-LC-L4-LCCS-Map-300m-P1Y-%d-v2.0.7cds.nc",
    plan_year
  ),
  sprintf(
    "C3S-LC-L4-LCCS-Map-300m-P1Y-%d-v2.1.1.nc",
    plan_year
  )
)
plan_path <- file.path(cci_dir, plan_name)

missing_files <- plan_path[!file.exists(plan_path)]
if (length(missing_files)) {
  stop(
    "Missing expected CCI NetCDF files:\n",
    paste(missing_files, collapse = "\n")
  )
}

out_tif <- file.path(out_dir, sprintf("ESACCI_frac_%d_0p05.tif", plan_year))


# Main loop

for (i in seq_along(plan_year)) {
  yr <- plan_year[i]
  f <- plan_path[i]
  ot <- out_tif[i]

  if (file.exists(ot)) {
    message("✓ Year ", yr, " already complete — skipping.")
    next
  }

  t0 <- Sys.time()
  message("→ [", yr, "] start")

  r <- rast(f, subds = "lccs_class")
  r <- terra::subst(r, nodata_vals, NA)

  m_stack <- rast(lapply(groups, function(cls) {
    # Classify the raster based on the current group
    classify(r, cbind(cls, 1), others = 0)
  }))

  frac <- resample(m_stack, tmpl, method = "average")
  names(frac) <- paste0("frac_", names(groups))

  w30 <- cfg$esa_cci$weights$cls30 %||% 0.75
  w40 <- cfg$esa_cci$weights$cls40 %||% 0.25

  fc <- frac$frac_cropland %||% 0
  fu <- frac$frac_urban %||% 0
  f30 <- frac$frac_cls30 %||% 0
  f40 <- frac$frac_cls40 %||% 0

  frac_fused <- clamp(fc + fu + w30 * f30 + w40 * f40, 0, 1)
  names(frac_fused) <- "frac_fused"

  frac_grass <- frac$frac_grassland
  names(frac_grass) <- "frac_grass"

  out <- c(frac_fused, frac_grass)
  names(out) <- c("frac_fused", "frac_grass")

  writeRaster(out,
    ot,
    overwrite = TRUE,
    gdal = gdal_opts,
    NAflag = -9999
  )

  rm(r, m_stack, frac, out)
  gc()

  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message("✓ [", yr, "] done (", dt, " s)")
}

message("Done.")
