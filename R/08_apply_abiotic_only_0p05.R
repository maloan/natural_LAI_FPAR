## =============================================================================
# 08_apply_abiotic_only_0p05.R — Apply abiotic (water/ice) drop mask to
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
VAR <- toupper(Sys.getenv("VAR", "FPAR"))
stopifnot(VAR %in% c("LAI", "FPAR"))

in_dir <- switch(VAR,
  LAI  = cfg$paths$georef_lai_0p05_dir,
  FPAR = cfg$paths$georef_fpar_0p05_dir
)


# Abiotic mask (FILE)
abi_mask_path <- file.path(
  cfg$paths$mask_abiotic_dir,
  "mask_abiotic_CCI_2007_tauW0p05_tauI0p05_0p05.tif"
)

if (is.null(abi_mask_path) || !file.exists(abi_mask_path)) {
  stop("Abiotic mask not found at cfg$paths$mask_abiotic_dir")
}
abi_mask <- rast(abi_mask_path)

# Output directory (DIR)
out_dir <- file.path(cfg$paths$out_root, "abiotic_only_0p05", VAR)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wopt <- wopt_f32(opts$SPEED_OVER_SIZE)$wopt
# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------
stopifnot(dir.exists(in_dir))
files <- sort(list.files(
  in_dir,
  pattern = paste0("^", VAR, "_\\d{6}_0p05\\.tif$"),
  full.names = TRUE
))
stopifnot(length(files) > 0L)

# ------------------------------------------------------------------------------
# Loop
# ------------------------------------------------------------------------------
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_%s_0p05_masked_abiotic.tif", VAR, ym))

  if (opts$SKIP_EXISTING && file.exists(out)) next

  r <- rast(f)
  if (!compareGeom(r, abi_mask, stopOnError = FALSE)) {
    r <- resample(r, abi_mask, method = "bilinear")
  }

  r <- mask(r, abi_mask, maskvalues = 1, updatevalue = NA)
  writeRaster(r, out, overwrite = TRUE, wopt = wopt)

  message("Wrote: ", out)
}

gc()
message(sprintf("Finished abiotic-only masking (0.05°): %s", VAR))
