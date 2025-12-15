## =============================================================================
# 08_abiotic_static_from_glc_0p05.R — Build static abiotic mask from GLC stack
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(stringr)
  library(glue)
  library(here)
})

# --- config & refs -------------------------------------------------------------
ROOT <- here()

source(here("R", "utils.R"))
source(here("R", "io.R"))
source(here("R", "geom.R"))
source(here("R", "viz.R"))
source(here("R", "options.R"))

cfg  <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

# --- references & paths --------------------------------------------------------
ref005  <- rast(cfg$grids$grid_005$ref_raster)
glc_dir <- cfg$paths$glc_dir

# --- thresholds & window -------------------------------------------------------
TAU_WATER <- as.numeric(Sys.getenv("TAU_WATER", "0.05"))
TAU_ICE   <- as.numeric(Sys.getenv("TAU_ICE", "0.05"))

Y1 <- as.integer(Sys.getenv("ABIOTIC_Y1", cfg$project$years$glc_start))
Y2 <- as.integer(Sys.getenv("ABIOTIC_Y2", cfg$project$years$glc_end))

# --- GLC class codes -----------------------------------------------------------
GLC        <- cfg$glc$classes
vals_water <- as.integer(unlist(GLC$water, use.names = FALSE))
vals_ice   <- as.integer(unlist(GLC$snow_ice, use.names = FALSE))
nodata_vals <- as.integer(unlist(GLC$nodata, use.names = FALSE))

# --- discover year files -------------------------------------------------------
files <- list.files(glc_dir, "\\.tif$", full.names = TRUE)
stopifnot(length(files) > 0)
yrs  <- as.integer(str_extract(basename(files), "(19|20)\\d{2}"))
keep <- which(!is.na(yrs) & yrs >= Y1 & yrs <= Y2)
files <- files[keep]
yrs   <- yrs[keep]
stopifnot(length(files) > 0)

# --- prepare accumulation ------------------------------------------------------
n_years <- 0L

out_dir <- here(cfg$paths$masks_root_dir, "mask_abiotic")
dir.create(out_dir, TRUE, showWarnings = FALSE)

out_tif <- here(
  cfg$paths$masks_root_dir,
  "mask_abiotic",
  sprintf(
    "mask_static_abiotic_GLC_%d-%d_tauW%s_tauI%s_tauB%s_0p05.tif",
    Y1,
    Y2,
    gsub("\\.", "p", sprintf("%.2f", as.numeric(TAU_WATER))),
    gsub("\\.", "p", sprintf("%.2f", as.numeric(TAU_ICE))),
  )
)

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
if (SKIP_EXISTING && file.exists(out_tif)) {
  message("✓ Abiotic GLC mask exists — skipping: ", out_tif)
  quit(save = "no")
}

# --- loop over years -----------------------------------------------------------
for (f in files) {
  r <- rast(f)

  if (is.na(crs(r))) {
    crs(r) <- crs(ref005)
  }

  if (!compareGeom(r, ref005, stopOnError = FALSE)) {
    r <- resample(r, ref005, method = "near")
  }

  if (length(nodata_vals)) {
    r[r %in% nodata_vals] <- NA
  }

  mW <- classify(r, cbind(vals_water, 1), others = 0)
  mI <- classify(r, cbind(vals_ice, 1), others = 0)

  n_years <- n_years + 1L
}

if (n_years == 0L) {
  stop("No GLC years found in window ", Y1, "-", Y2)
}

# --- averages ------------------------------------------------------------------
pW <- mW / n_years
pI <- mI / n_years

# --- threshold + combine -------------------------------------------------------
abiotic  <- (pW >= TAU_WATER) | (pI >= TAU_ICE)
abi_mask <- ifel(abiotic, 1L, 0L)
names(abi_mask) <- "abiotic_drop"

# --- write output --------------------------------------------------------------
writeRaster(
  abi_mask,
  out_tif,
  overwrite = TRUE,
  gdal = gdal_wopt("LOG1S")$gdal,
  NAflag = 255
)

cat(
  glue(
    "
Wrote abiotic overlay: {out_tif}
Window: {Y1}-{Y2}  (n={n_years} years)
Thresholds — water={TAU_WATER}, ice={TAU_ICE}
"
  )
)
gc()
