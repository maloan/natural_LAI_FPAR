## =============================================================================
# 02_apply_abiotic_only_0p05.R — Apply abiotic (water/ice) drop mask to
# georeferenced monthly LAI/FPAR at 0.05°
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------
var <- toupper(Sys.getenv("var", "FPAR"))
stopifnot(var %in% c("LAI", "FPAR"))

var_lower <- tolower(var)
in_dir <- cfg$paths[[sprintf("georef_%s_0p05_dir", var_lower)]]

# non-vegetated mask
abi_mask_path <- file.path(
  cfg$paths$mask_abiotic_dir,
  "mask_abiotic_CCI_2007_tauW0p05_tauI0p05_0p05.tif"
)
abi_mask <- rast(abi_mask_path)

# Output directory
out_dir <- file.path(cfg$paths$out_root, "abiotic_only_0p05", var)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wopt <- wopt_f32(opts$speed_over_size)
# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------
files <- sort(list.files(
  in_dir,
  pattern = paste0("^", var, "_\\d{6}_0p05\\.tif$"),
  full.names = TRUE
))

# ------------------------------------------------------------------------------
# Loop
# ------------------------------------------------------------------------------
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_%s_0p05_masked_abiotic.tif", var, ym))

  if (opts$skip_existing && file.exists(out)) {
    next
  }

  r <- rast(f)
  if (!compareGeom(r, abi_mask, stopOnError = FALSE)) {
    r <- resample(r, abi_mask, method = "bilinear")
  }

  r <- mask(r, abi_mask, maskvalues = 1, updatevalue = NA)
  writeRaster(r, out, overwrite = TRUE, wopt = wopt)

  message("Wrote: ", out)
}

gc()
message(sprintf("Finished abiotic-only masking (0.05°): %s", var))
