## =============================================================================
# 01_georef_0p05.R — Georeference monthly LAI/FPAR to 0.05° grid
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "geom.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------
VAR <- toupper(Sys.getenv("VAR", "FPAR"))
stopifnot(VAR %in% c("LAI", "FPAR"))

REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))
wopt <- wopt_f32(opts$SPEED_OVER_SIZE)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
Vcfg <- cfg$variables[[tolower(VAR)]]

paths <- switch(VAR,
  LAI = list(
    out_dir = cfg$paths$georef_lai_0p05_dir, in_dir = cfg$paths$lai_nc_dir
  ),
  FPAR = list(
    out_dir = cfg$paths$georef_fpar_0p05_dir, in_dir = cfg$paths$fpar_nc_dir
  )
)

out_georef <- paths$out_dir
out_quick <- file.path(out_georef, "quicklooks")
nc_dir <- paths$in_dir

dir.create(out_georef, recursive = TRUE, showWarnings = FALSE)
dir.create(out_quick, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------
stopifnot(dir.exists(nc_dir))
files <- sort(list.files(nc_dir, pattern = "\\.nc(4)?$", full.names = TRUE))
stopifnot(length(files) > 0L)

# ------------------------------------------------------------------------------
# Loop
# ------------------------------------------------------------------------------
for (f in files) {
  ym <- extract_ym_from_filename(f)

  out_tif <- file.path(out_georef, sprintf("%s_%s_0p05.tif", VAR, ym))
  out_png <- file.path(out_quick, sprintf("%s_%s_0p05.png", VAR, ym))

  if (!opts$FORCE && opts$SKIP_EXISTING && file.exists(out_tif)) next

  r <- nc_month_to_raster(
    nc_file     = f,
    ym          = ym,
    VAR         = VAR,
    Vcfg        = Vcfg,
    ref         = ref005,
    method      = "bilinear",
    strict_time = TRUE
  )

  writeRaster(r, out_tif, overwrite = TRUE)

  if (REMAKE_QL || !file.exists(out_png)) {
    write_quicklook_raster(
      r, out_png,
      title = sprintf("%s %s (0.05°)", VAR, ym)
    )
  }
}

gc()
message("Done georeferencing (0.05°).")
