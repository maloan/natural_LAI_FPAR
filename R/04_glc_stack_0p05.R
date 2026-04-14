## =============================================================================
# 04_glc_stack_0p05.R — Build annual GLC_FCS30D categorical yearstack (0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "paths.R"))
source(here("R", "helpers", "files.R"))
source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "plotting.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.9)

ref005 <- rast(cfg$grids$grid_005$ref_raster)

glc_dir <- cfg$paths$glc_dir
out_dir <- cfg$paths$glc_out_dir
stack_out <- file.path(out_dir, "glc_cat_yearstack_0p05.tif")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

skip_existing <- as_bool(Sys.getenv("skip_existing"), default = TRUE)
overwrite <- as_bool(Sys.getenv("overwrite"), default = FALSE)

cropland_vals <- as.integer(unlist(cfg$glc$classes$cropland))
urban_vals <- as.integer(unlist(cfg$glc$classes$urban))
nodata_vals <- as.integer(unlist(cfg$glc$classes$nodata))

years_for_ql <- intersect(c(1990, 2000, 2010, 2020), cfg$glc$years)

# --- discover files (base R) --------------------------------------------------
all_files <- list.files(glc_dir, pattern = "\\.tif$", full.names = TRUE)
if (!length(all_files)) {
  stop("No GLC GeoTIFFs found in: ", glc_dir)
}

# Extract years from filenames (first 4 characters are year YYYY)
basenames <- basename(all_files)
years <- as.integer(substr(basenames, 1, 4))

# Filter to valid years and sort
valid_mask <- years %in% cfg$glc$years
ord <- order(years[valid_mask])

paths <- all_files[valid_mask][ord]
years <- years[valid_mask][ord]

stopifnot(length(paths) > 0)

message(
  "Found ", length(paths),
  " GLC rasters; span [", min(years), "..", max(years), "]."
)

# --- rebuild stack ------------------------------------------------------------
bands <- vector("list", length(paths))

for (i in seq_along(paths)) {
  yr <- years[i]
  f <- paths[i]

  message("→ Processing ", basename(f))

  r <- rast(f)[[1]]
  if (is.na(crs(r))) {
    crs(r) <- crs(ref005)
  }

  if (!compareGeom(r, ref005, stopOnError = FALSE)) {
    r <- resample(r, ref005, method = "near")
  }

  if (length(nodata_vals)) {
    r <- terra::subst(r, nodata_vals, NA)
  }

  names(r) <- sprintf("Y%04d", yr)
  bands[[i]] <- r

  gc()
}

stack <- rast(bands)
writeRaster(stack, stack_out, overwrite = TRUE)
gc()

message("Wrote GLC yearstack: ", stack_out)
