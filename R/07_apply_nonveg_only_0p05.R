## =============================================================================
# Step 07 — Apply non-vegetated (water/ice) drop mask to
# georeferenced monthly LAI/FPAR at 0.05°
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "paths.R"))
source(here("R", "helpers", "files.R"))
source(here("R", "helpers", "netcdf.R"))
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
nonveg_mask_path <- file.path(
  cfg$paths$mask_nonvegetated_dir,
  "mask_nonvegetated_CCI_2007_tauW0p05_tauI0p05_0p05.tif"
)
nonveg_mask <- rast(nonveg_mask_path)

# Output directory
out_dir <- file.path(here("output"), "nonvegetated_only_0p05", var)
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
  out <- file.path(out_dir, sprintf("%s_%s_0p05_masked_nonvegetated.tif", var, ym))

  if (opts$skip_existing && file.exists(out)) {
    next
  }

  r <- rast(f)
  r <- align_to_template(r, nonveg_mask, method = "bilinear")

  r <- mask(r, nonveg_mask, maskvalues = 1, updatevalue = NA)
  writeRaster(r, out, overwrite = TRUE, wopt = wopt)

  message("Wrote: ", out)
}

gc()
message(sprintf("Finished non-vegetated-only masking (0.05°): %s", var))
