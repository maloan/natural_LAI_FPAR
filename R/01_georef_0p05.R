# =============================================================================
# 01_georef_0p05.R — Georeference monthly LAI/FPAR to 0.05° grid
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

# Config
var <- toupper(Sys.getenv("var", "FPAR"))
stopifnot(var %in% c("LAI", "FPAR"))

var_lower <- tolower(var)
wopt <- wopt_f32(opts$speed_over_size)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
vcfg <- cfg$variables[[var_lower]]

# Dynamically construct paths using config keys
out_georef <- cfg$paths[[sprintf("georef_%s_0p05_dir", var_lower)]]
nc_dir <- cfg$paths[[sprintf("%s_nc_dir", var_lower)]]

lapply(c(out_georef),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

# Inputs
stopifnot(dir.exists(nc_dir))
files <- sort(list.files(nc_dir, pattern = "\\.nc(4)?$", full.names = TRUE))
stopifnot(length(files) > 0L)

# Loop
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out_tif <- file.path(out_georef, sprintf("%s_%s_0p05.tif", var, ym))
  if (!opts$force && opts$skip_existing && file.exists(out_tif)) {
    next
  }

  r <- nc_month_to_raster(
    nc_file     = f,
    ym          = ym,
    var         = var,
    vcfg        = vcfg,
    ref         = ref005,
    method      = "bilinear",
    strict_time = TRUE
  )
  writeRaster(r, out_tif, overwrite = TRUE)
}

gc()
