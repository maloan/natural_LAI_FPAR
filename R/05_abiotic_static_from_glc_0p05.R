## =============================================================================
# 05_abiotic_static_from_glc_0p05.R — Build static abiotic mask from GLC stack
#
# Purpose
#   Create a static “abiotic” drop mask (1=drop, 0=keep) using multi-year
#   categorical GLC_FCS30D maps aggregated to 0.05° resolution, marking pixels
#   persistently dominated by water, snow/ice, or bare surfaces.
#
# Inputs
#   - Annual GLC_FCS30D GeoTIFFs (0.05°): cfg$paths$glc_dir
#   - Reference grid (0.05°): cfg$grids$grid_005$ref_raster
#   - Class codes: cfg$glc$classes (water, snow_ice, bare, nodata)
#
# Outputs
#   - mask_static_abiotic_GLC_{Y1–Y2}_tauW{τW}_tauI{τI}_tauB{τB}_0p05.tif
#     (Byte; 1=drop, 0=keep, NA=255)
#
# Environment variables
#   TAU_WATER (numeric, default 0.05)
#   TAU_ICE   (numeric, default 0.05)
#   TAU_BARE  (numeric, default 0.30)
#   ABIOTIC_Y1, ABIOTIC_Y2 (integer; year window)
#   SKIP_EXISTING (logical, default TRUE)
#
# Dependencies
#   Packages: terra, stringr, glue
#   Sourced helpers: utils.R, io.R, geom.R, viz.R, options.R
#
# Processing overview
#   1) Identify yearly GLC rasters within the configured time window.
#   2) Build binary masks (0/1) for water, ice, and bare classes.
#   3) Average across years to obtain multi-year class fractions.
#   4) Drop pixels where pW≥τW ∨ pI≥τI ∨ pB≥τB.
#   5) Write final static abiotic mask (Byte; 1=drop, 0=keep).
## =============================================================================


suppressPackageStartupMessages({
  library(terra)
  library(stringr)
  library(glue)
})

# --- config & refs -------------------------------------------------------------
# Root dir
ROOT <- normalizePath(path.expand(
  Sys.getenv("SNU_LAI_FPAR_ROOT", "~/GitHub/natural_LAI_FPAR")
),
winslash = "/",
mustWork = FALSE)
# Other source files
source(file.path(ROOT, "R", "utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))
source(file.path(ROOT, "R", "viz.R"))
source(file.path(ROOT, "R", "options.R"))
cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)

# --- references & paths --------------------------------------------------------
ref005  <- rast(cfg$grids$grid_005$ref_raster)
glc_dir <- cfg$paths$glc_dir   # 0.05° annual categorical GeoTIFFs

# --- thresholds & window -------------------------------------------------------
TAU_WATER <- as.numeric(Sys.getenv("TAU_WATER", "0.05"))
TAU_ICE   <- as.numeric(Sys.getenv("TAU_ICE", "0.05"))
TAU_BARE  <- as.numeric(Sys.getenv("TAU_BARE", "0.30"))

Y1 <- as.integer(Sys.getenv("ABIOTIC_Y1", cfg$project$years$glc_start))
Y2 <- as.integer(Sys.getenv("ABIOTIC_Y2", cfg$project$years$glc_end))

tok <- function(x)
  gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))  # for filenames

# --- GLC class codes -----------------------------------------------------------
.as_int    <- function(x)
  as.integer(unlist(x, use.names = FALSE))
GLC        <- cfg$glc$classes
vals_water <- .as_int(GLC$water)
vals_ice   <- .as_int(GLC$snow_ice)
vals_bare  <- .as_int(GLC$bare)
nodata_vals <- .as_int(GLC$nodata)

# --- discover year files -------------------------------------------------------
files <- list.files(glc_dir, "\\.tif$", full.names = TRUE)
stopifnot(length(files) > 0)

yrs <- as.integer(str_extract(basename(files), "(19|20)\\d{2}"))
keep <- which(!is.na(yrs) & yrs >= Y1 & yrs <= Y2)
files <- files[keep]
yrs <- yrs[keep]
stopifnot(length(files) > 0)

# --- accumulate per-class means over years ------------------------------------
sumW <- sumI <- sumB <- NULL
n_years <- 0L

bin <- function(rr, codes)
  classify(rr, cbind(codes, 1), others = 0)

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
if (SKIP_EXISTING && file.exists(out_tif)) {
  message("✓ Abiotic GLC mask exists — skipping: ", out_tif)
  quit(save = "no")
}


for (f in files) {
  r <- rast(f)
  if (is.na(crs(r))) {
    crs(r) <- crs(ref005)
  }
  if (!compareGeom(r, ref005, stopOnError = FALSE)) {
    r <- resample(r, ref005, method = "near")  # categorical
  }
  if (length(nodata_vals))
  {
    r[r %in% nodata_vals] <- NA
  }

  # 0/1 masks on the 0.05° grid
  mW <- bin(r, vals_water)
  mI <- bin(r, vals_ice)
  mB <- bin(r, vals_bare)

  # accumulate; mean across years will be fractions in [0,1]
  sumW <- if (is.null(sumW)) {
    mW
  } else {
    sumW + mW
  }
  sumI <- if (is.null(sumI)) {
    mI
  } else {
    sumI + mI
  }
  sumB <- if (is.null(sumB)) {
    mB
  } else {
    sumB + mB
  }
  n_years <- n_years + 1L
}

if (n_years == 0L) {
  stop("No GLC years found in window ", Y1, "-", Y2)
}

# --- average across years ------------------------------------------------------
pW <- sumW / n_years
pI <- sumI / n_years
pB <- sumB / n_years

# --- threshold + combine (OR) --------------------------------------------------
abiotic  <- (pW >= TAU_WATER) | (pI >= TAU_ICE) | (pB >= TAU_BARE)
abi_mask <- ifel(abiotic, 1L, 0L)  # 1=drop, 0=keep
names(abi_mask) <- "abiotic_drop"

# --- write --------------------------------------------------------------------
out_dir <- file.path(cfg$paths$masks_root_dir, "mask_abiotic")
dir.create(out_dir, TRUE, showWarnings = FALSE)
out_tif <- file.path(
  out_dir,
  sprintf(
    "mask_static_abiotic_GLC_%d-%d_tauW%s_tauI%s_tauB%s_0p05.tif",
    Y1,
    Y2,
    tok(TAU_WATER),
    tok(TAU_ICE),
    tok(TAU_BARE)
  )
)

writeRaster(
  abi_mask,
  out_tif,
  overwrite = TRUE,
  gdal = gdal_wopt("LOG1S")$gdal,
  NAflag = 255
)
gc()
cat(
  glue(
    "
Wrote abiotic overlay: {out_tif} window: {Y1}-{Y2} (n={n_years} years) thresholds: water={TAU_WATER}, ice={TAU_ICE}, bare={TAU_BARE}
"
  )
)
