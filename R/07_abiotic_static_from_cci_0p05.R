## =============================================================================
# 07_abiotic_static_from_cci.R — Build abiotic mask for ONE year (0.05°)
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

terraOptions(progress = 1, memfrac = 0.6)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
cci_dir <- cfg$paths$cci_dir

tau_water <- as.numeric(Sys.getenv("tau_water", "0.05"))
tau_ice <- as.numeric(Sys.getenv("tau_ice", "0.05"))
skip_existing <- as_bool(Sys.getenv("skip_existing"), default = FALSE)
year <- as.integer(Sys.getenv("year", "2007"))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_abiotic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tau_w_tok <- gsub("\\.", "p", sprintf("%.2f", tau_water))
tau_i_tok <- gsub("\\.", "p", sprintf("%.2f", tau_ice))

out_tif <- file.path(out_dir, sprintf(
  "mask_abiotic_CCI_%d_tauW%s_tauI%s_0p05.tif", year, tau_w_tok, tau_i_tok
))

if (skip_existing && file.exists(out_tif)) {
  message("✓ Abiotic mask already exists — skipping: ", out_tif)
  return(invisible(NULL))
}

esa_cci <- cfg$esa_cci$classes
vals_water <- as.integer(unlist(esa_cci$water))
vals_ice <- as.integer(unlist(esa_cci$snow_ice))
nodata_vals <- unique(c(as.integer(unlist(esa_cci$nodata)), 255L))

# --- input for year ---
files <- list.files(cci_dir, pattern = "\\.tif$", full.names = TRUE)
basenames <- basename(files)

# Extract years from filenames
yrs <- as.integer(substr(basenames, 1, 4))
cand <- files[!is.na(yrs) & yrs == year]

# Rank by source (C3S preferred = 2, others = 1)
cand_rank <- ifelse(grepl("^C3S", basename(cand)), 2L, 1L)
in_file <- cand[which.max(cand_rank)]

message("→ Processing ", basename(in_file), " (year=", year, ")")

r <- rast(in_file)
if (is.na(crs(r))) {
  crs(r) <- crs(ref005)
}

r <- terra::subst(r, nodata_vals, NA)

# Water and ice fractions at 0.05° (single year)
pW <- resample(classify(r, cbind(vals_water, 1), others = 0), ref005, method = "average")
pI <- resample(classify(r, cbind(vals_ice, 1), others = 0), ref005, method = "average")

# Create drop masks
water_drop <- ifel(pW >= tau_water, 1L, 0L)
ice_drop <- ifel(pI >= tau_ice, 1L, 0L)
both_drop <- ifel(water_drop & ice_drop, 1L, 0L)
abi_mask_combined <- ifel(water_drop | ice_drop, 1L, 0L)

# Write component rasters
writeRaster(c(water_drop, ice_drop, both_drop, abi_mask_combined),
  file.path(out_dir, sprintf("abiotic_components_%d_0p05.tif", year)),
  overwrite = TRUE, wopt = wopt_byte(FALSE, na = 255L)
)

# Write final combined mask
names(abi_mask_combined) <- "abiotic_drop"
writeRaster(
  abi_mask_combined,
  out_tif,
  overwrite = TRUE,
  wopt = wopt_byte(FALSE, na = 255L)
)

gc()
