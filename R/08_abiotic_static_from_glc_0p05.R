## =============================================================================
# 08_abiotic_static_from_glc_0p05.R — Build static abiotic mask from GLC years
#   Semantics: 1 = drop (abiotic), 0 = keep, 255 = NA
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.25)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
glc_dir <- cfg$paths$glc_dir

TAU_WATER <- as.numeric(Sys.getenv("TAU_WATER", "0.05"))
TAU_ICE <- as.numeric(Sys.getenv("TAU_ICE", "0.05"))

Y1 <- as.integer(Sys.getenv("ABIOTIC_Y1", cfg$project$years$glc_start))
Y2 <- as.integer(Sys.getenv("ABIOTIC_Y2", cfg$project$years$glc_end))

GLC <- cfg$glc$classes
vals_water <- as.integer(unlist(GLC$water, use.names = FALSE))
vals_ice <- as.integer(unlist(GLC$snow_ice, use.names = FALSE))
nodata_vals <- as.integer(unlist(GLC$nodata, use.names = FALSE))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_abiotic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tok <- function(x) gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
out_tif <- file.path(out_dir, sprintf(
  "mask_static_abiotic_GLC_%d-%d_tauW%s_tauI%s_0p05.tif",
  Y1, Y2, tok(TAU_WATER), tok(TAU_ICE)
))

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
if (SKIP_EXISTING && file.exists(out_tif)) {
  message("✓ Abiotic GLC mask exists — skipping: ", out_tif)
  quit(save = "no")
}

# --- discover yearly rasters ---------------------------------------------------
files <- list.files(glc_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(files) > 0L)

get_year <- function(p) {
  m <- regexpr("(19|20)\\d{2}", basename(p), perl = TRUE)
  if (m[1] < 0) {
    return(NA_integer_)
  }
  as.integer(substr(basename(p), m[1], m[1] + attr(m, "match.length") - 1))
}
yrs <- vapply(files, get_year, integer(1))

keep <- which(!is.na(yrs) & yrs >= Y1 & yrs <= Y2)
files <- files[keep]
yrs <- yrs[keep]
stopifnot(length(files) > 0L)

# --- accumulate counts across years -------------------------------------------
sumW <- NULL
sumI <- NULL
n_years <- 0L

for (f in files) {
  r <- rast(f)[[1]]

  if (is.na(crs(r))) crs(r) <- crs(ref005)
  if (!compareGeom(r, ref005, stopOnError = FALSE)) {
    r <- resample(r, ref005, method = "near")
  }
  if (length(nodata_vals)) r[r %in% nodata_vals] <- NA

  mW <- classify(r, cbind(vals_water, 1), others = 0)
  mI <- classify(r, cbind(vals_ice, 1), others = 0)

  if (is.null(sumW)) {
    sumW <- mW
    sumI <- mI
  } else {
    sumW <- sumW + mW
    sumI <- sumI + mI
  }

  n_years <- n_years + 1L
  gc()
}

stopifnot(n_years > 0L)

pW <- sumW / n_years
pI <- sumI / n_years

abi_mask <- ifel((pW >= TAU_WATER) | (pI >= TAU_ICE), 1L, 0L)
names(abi_mask) <- "abiotic_drop"

writeRaster(
  abi_mask,
  out_tif,
  overwrite = TRUE,
  gdal = gdal_wopt("LOG1S")$gdal,
  NAflag = 255
)

cat(sprintf(
  "Wrote: %s\nWindow: %d-%d (n=%d)\nThresholds: water=%.2f, ice=%.2f\n",
  out_tif, Y1, Y2, n_years, TAU_WATER, TAU_ICE
))
