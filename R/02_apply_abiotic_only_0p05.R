## =============================================================================
# 02_apply_abiotic_only_0p05.R — Remove water and ice from the georef data -
# Manually change between LAI and FPAR
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(stringr)
  library(glue)
  library(rlist)
})

# georef LAI directory (0.05°)
in_dir  <- "~/GitHub/natural_LAI_FPAR/data/georef/georef_lai_0p05"

# abiotic mask (0.05°), values 0 = keep, 1 = drop
abi_mask_path <- "~/GitHub/natural_LAI_FPAR/output/tau_0.1/masks/mask_abiotic/mask_static_abiotic_CCI_tauW0p05_tauI0p05_0p05.tif"

# output directory
out_dir <- "~/GitHub/natural_LAI_FPAR/output/georef_biotic/georef_biotic_LAI_0p05"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------
## Load reference mask
## ------------------------------------------------------------
abi_mask <- rast(abi_mask_path)
ref <- abi_mask

## ------------------------------------------------------------
## Find all georef LAI files
## ------------------------------------------------------------
files <- list.files(in_dir, "^LAI_\\d{6}_0p05\\.tif$", full.names = TRUE)
files <- list.reverse(files)
if (!length(files)) {
  stop("No LAI files found.")
}

## ------------------------------------------------------------
## Process each file
## ------------------------------------------------------------
for (f in files) {
  ym <- str_extract(basename(f), "\\d{6}")
  out <- file.path(out_dir, glue("LAI_{ym}_0p05_masked_abiotic.tif"))
  r <- rast(f)

  # regrid if needed
  if (!compareGeom(r, ref, stopOnError = FALSE)) {
    r <- resample(r, ref, method = "bilinear")
  }

  # mask: abi_mask==1 → NA
  r_masked <- mask(r, abi_mask, maskvalues = 1, updatevalue = NA)

  writeRaster(
    r_masked,
    out,
    overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"),
    NAflag = -9999
  )

  message("Wrote: ", out)
}

message("=== Finished abiotic-only masking of LAI 0.05° ===")
