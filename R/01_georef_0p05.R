# =============================================================================
# 01_georef_0p05.R — Georeference monthly LAI/FPAR to 0.05° grid
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.25)

# Config
var <- toupper(Sys.getenv("VAR"))
stopifnot(var %in% c("LAI", "FPAR"))

var_lower <- tolower(var)
wopt <- wopt_f32(FALSE)

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
months <- format(
  seq(
    as.Date(sprintf("%d-01-01", cfg$project$years$lai_start)),
    as.Date(sprintf("%d-12-01", cfg$project$years$lai_end)),
    by = "month"
  ),
  "%Y%m"
)
input_months <- vapply(files, extract_ym_from_filename, character(1),
                       USE.NAMES = FALSE)
files <- files[order(input_months)]

dependencies <- c(
  cfg$grids$grid_005$ref_raster,
  here("config", sprintf("config_%s.yml", cfg$project$run_tag)),
  here("R", "01_georef_0p05.R"),
  here("R", "helpers", "netcdf.R"),
  here("R", "helpers", "io.R")
)
stopifnot(all(file.exists(dependencies)))
source_mtime <- max(file.mtime(dependencies))

# Loop
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out_tif <- file.path(out_georef, sprintf("%s_%s_0p05.tif", var, ym))

  r <- nc_month_to_raster(
    nc_file     = f,
    ym          = ym,
    var         = var,
    vcfg        = vcfg,
    ref         = ref005,
    method      = "bilinear",
    strict_time = TRUE
  )
  tmp <- paste0(out_tif, ".tmp.tif")
  writeRaster(r, tmp, overwrite = TRUE, wopt = wopt)
  if (!file.rename(tmp, out_tif)) {
    stop("Could not replace ", out_tif)
  }
}

gc()
