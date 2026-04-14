## =============================================================================
# 11_agg_0p5.R — Area-weighted aggregation of masked LAI/FPAR (0.05° → 0.5°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

# --- repo helpers --------------------------------------------------------------
source(here("R", "helpers", "paths.R"))
source(here("R", "helpers", "files.R"))
source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

# ------------------------------------------------------------------------------
# Env
# ------------------------------------------------------------------------------
skip_existing <- as_bool(Sys.getenv("skip_existing"), default = FALSE)
overwrite <- as_bool(Sys.getenv("overwrite"), default = FALSE)

var <- toupper(Sys.getenv("var", "FPAR")) # LAI|FPAR
mask <- toupper(Sys.getenv("mask", "CCI")) # CCI|GLC

# ------------------------------------------------------------------------------
# Refs and weights
# ------------------------------------------------------------------------------
ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref050 <- aggregate(ref005, fact = 10, fun = "mean", na.rm = TRUE) # template
area005 <- rast(cfg$grids$grid_005$area_raster)

# ------------------------------------------------------------------------------
# dirs
# ------------------------------------------------------------------------------
key_in <- sprintf("masked_%s_%s_0p05_dir", tolower(var), tolower(mask))
key_out <- sprintf("masked_%s_%s_0p5_dir", tolower(var), tolower(mask))

in_dir <- cfg$paths[[key_in]]
out_dir <- cfg$paths[[key_out]]

stopifnot(is.character(in_dir), length(in_dir) == 1, nzchar(in_dir), dir.exists(in_dir))
stopifnot(is.character(out_dir), length(out_dir) == 1, nzchar(out_dir))
patt <- sprintf("^%s_\\d{6}_0p05_masked\\.tif$", var)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# inputs
# ------------------------------------------------------------------------------
stopifnot(dir.exists(in_dir))
files <- sort(list.files(in_dir, pattern = patt, full.names = TRUE))
stopifnot(length(files) > 0L)

# --- loop ---
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_masked_%s_0p5.tif", var, ym))
  do_write <- overwrite || !file.exists(out) || !skip_existing
  if (!do_write) next

  r <- rast(f)
  if (!compareGeom(r, ref005, stopOnError = FALSE)) {
    r <- resample(r, ref005, method = "bilinear")
  }

  num <- aggregate(r * area005, fact = 10, fun = "sum", na.rm = TRUE)
  den <- aggregate((!is.na(r)) * area005, fact = 10, fun = "sum", na.rm = TRUE)
  r050 <- ifel(den == 0, NA, num / den)

  if (!compareGeom(r050, ref050, stopOnError = FALSE)) {
    r050 <- resample(r050, ref050, method = "near")
  }

  wopt <- wopt_f32(opts$speed_over_size)
  writeRaster(r050, out, overwrite = TRUE, wopt = wopt)

  gc()
}

message("Aggregated monthly ", var, " written to: ", out_dir)
