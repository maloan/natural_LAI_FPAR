#!/usr/bin/env Rscript
## =============================================================================
# 04_luh_pasture_overlap_025.R — Build LUH pasture-overlap mask (0.25° + 0.05°)
#
# Purpose
#   Identify regions where LUH pasture fraction overlaps substantially with
#   satellite-derived grass cover. Outputs a binary “drop” mask:
#     1 = drop (likely managed grass), 0 = keep, NA = nodata (written as 255).
#
# Inputs
#   - LUH2 states.nc (cfg$luh2$states_nc; subds = cfg$luh2$variables$pasture)
#   - Grass fraction at 0.05° from GRASS_SOURCE = {CCI | GLC_FRAC | GLC_TEMP}
#   - Reference grids: cfg$grids$grid_005$ref_raster, cfg$grids$grid_025$ref_raster
#   - Cell area raster (0.05°): cfg$grids$grid_005$area_raster
#
# Outputs
#   - {cfg$paths$masks_root_dir}/mask_luh_overlap/mask_luh_overlap_{SRC}_Gmin{G}_Pmin{P}_alpha{A}_{Y0}-{Y1}_0p25.tif
#   - {same dir}/mask_luh_overlap_{SRC}_Gmin{G}_Pmin{P}_alpha{A}_{Y0}-{Y1}_0p05_rep.tif
#   - Quicklook PNG (optional)
#
# Environment variables
#   GRASS_SOURCE (default "GLC_FRAC")   # {CCI | GLC_FRAC | GLC_TEMP}
#   ALPHA        (default 0.50)         # threshold on pasture/grass
#   LUH_AVG_START, LUH_AVG_END (defaults from cfg$project$years$cci_start/end)
#   REMAKE_QL    (default TRUE)
#
# Processing overview
#   1) Mean grass at 0.05° over [Y0..Y1] from source.
#   2) Area-weighted aggregation to 0.25°.
#   3) Mean LUH pasture at 0.25° over [Y0..Y1].
#   4) Mask rule: drop if (pasture/grass)≥ALPHA.
#   5) Write 0.25° mask and 0.05° replica; optionally quicklooks.
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(scico)
  library(yaml)
  library(glue)
  library(stringr)
})

## --- config & helpers ---------------------------------------------------------
ROOT <- normalizePath(path.expand(
  Sys.getenv("SNU_LAI_FPAR_ROOT", "~/GitHub/natural_LAI_FPAR")
), winslash = "/", mustWork = FALSE)

source(file.path(ROOT, "R", "utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))
source(file.path(ROOT, "R", "viz.R"))
source(file.path(ROOT, "R", "options.R"))

cfg  <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.30)

ref005  <- rast(cfg$grids$grid_005$ref_raster)
ref025  <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)

# Fallback for GDAL options helper if not present
safe_gdal <- function(profile = "LOG1S") {
  if (exists("gdal_wopt", mode = "function")) {
    return(gdal_wopt(profile)$gdal)
  } else {
    # LOG1S: logical/byte with compression; else default to Byte + LZW
    if (identical(profile, "LOG1S")) {
      return(c("COMPRESS=LZW", "PREDICTOR=2", "ZLEVEL=6", "TILED=YES"))
    } else {
      return(character(0))
    }
  }
}

## --- params (env) -------------------------------------------------------------
GRASS_SOURCE <- toupper(Sys.getenv("GRASS_SOURCE", "CCI"))  # CCI|GLC_TEMP
REMAKE_QL    <- as.logical(Sys.getenv("REMAKE_QL", "TRUE"))
ALPHA        <- as.numeric(Sys.getenv("ALPHA",  "0.50"))
Y0 <- env_get_int("LUH_AVG_START", cfg$project$years$cci_start)
Y1 <- env_get_int("LUH_AVG_END",   cfg$project$years$cci_end)

stopifnot(is.finite(ALPHA), is.finite(Y0), is.finite(Y1))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_luh_overlap")
ql_dir  <- file.path(out_dir, "quicklooks")
dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)

# Output basename
fname <- glue("mask_luh_overlap_{GRASS_SOURCE}_",
              "alpha{tok(ALPHA)}_{Y0}-{Y1}_0p25.tif")
pal_green <- hcl.colors(64, "Greens", rev = TRUE)

## --- 1) Grass fraction @0.05° over [Y0..Y1] ----------------------------------
grass_005 <- switch(
  GRASS_SOURCE,
  "CCI" = {
    frac_dir <- cfg$paths$cci_out_dir
    f <- list.files(frac_dir, "^ESACCI_frac_\\d{4}_0p05\\.tif$", full.names = TRUE)
    if (!length(f)) stop("No CCI fraction files in: ", frac_dir)
    yrs  <- as.integer(str_extract(basename(f), "\\d{4}"))
    keep <- which(yrs >= Y0 & yrs <= Y1)
    if (!length(keep)) stop("No CCI fraction years in window ", Y0, "-", Y1)
    stk <- rast(lapply(f[keep], \(x) rast(x)[["frac_grass"]]))
    mean(stk, na.rm = TRUE)
  },
  "GLC_TEMP" = {
    p <- file.path(cfg$paths$glc_out_dir, "glc_cat_yearstack_0p05.tif")
    if (!file.exists(p)) stop("GLC yearstack not found: ", p)
    s <- rast(p)
    y <- suppressWarnings(as.integer(substr(names(s), 2, 5)))
    keep <- which(y >= Y0 & y <= Y1)
    if (!length(keep)) stop("No GLC years in window ", Y0, "-", Y1)
    grass_vals <- as.integer(unlist(cfg$glc$classes$grassland, use.names = FALSE))
    is_grass <- classify(s[[keep]], cbind(grass_vals, 1), others = 0)
    app(is_grass, mean, na.rm = TRUE)
  },
  stop("Unknown GRASS_SOURCE='", GRASS_SOURCE, "'. Use CCI | GLC_FRAC | GLC_TEMP.")
)

if (!isTRUE(compareGeom(grass_005, ref005, stopOnError = FALSE))) {
  grass_005 <- resample(grass_005, ref005, method = "bilinear")
}
names(grass_005) <- "grass"

## --- 2) Aggregate grass 0.05° → 0.25° (area-weighted mean) -------------------
# Factor 5 assumed for 0.05 → 0.25; using area-weighted sum/valid-area sum.
num <- aggregate(grass_005 * area005, fact = 5, fun = \(x) sum(x, na.rm = TRUE))
den <- aggregate((!is.na(grass_005)) * area005, fact = 5, fun = \(x) sum(x, na.rm = TRUE))
grass_025 <- ifel(den == 0, NA, num / den)
if (!isTRUE(compareGeom(grass_025, ref025, stopOnError = FALSE))) {
  grass_025 <- resample(grass_025, ref025, method = "near")
}
names(grass_025) <- "grass_025"

## --- 3) LUH pasture @0.25° (time mean) ---------------------------------------
luh_nc <- cfg$luh2$states_nc
if (!file.exists(luh_nc)) stop("LUH file not found: ", luh_nc)
v_pas <- cfg$luh2$variables$pasture

# Read subdataset; silence CF bounds warnings from sds().
suppressWarnings({
  sds_names <- try(names(sds(luh_nc)), silent = TRUE)
})
if (inherits(sds_names, "try-error") || !(v_pas %in% sds_names)) {
  # fallback to NETCDF: path syntax
  pas <- rast(sprintf("NETCDF:%s:%s", luh_nc, v_pas))
} else {
  pas <- rast(luh_nc, subds = v_pas)
}

ty <- suppressWarnings(as.integer(time(pas)))
if (all(is.na(ty))) stop("LUH time axis missing/NA for '", v_pas, "'.")
keep <- which(ty >= Y0 & ty <= Y1)
if (!length(keep)) stop("No LUH timesteps in window ", Y0, "-", Y1)

pasture_025 <- mean(pas[[keep]], na.rm = TRUE)
if (!isTRUE(compareGeom(pasture_025, ref025, stopOnError = FALSE))) {
  pasture_025 <- resample(pasture_025, ref025, method = "bilinear")
}
pasture_025 <- clamp(pasture_025, 0, 1)
names(pasture_025) <- "pasture"

## --- 4) Decision rule at 0.25° -----------------------------------------------
eps   <- 1e-9
denom <- ifel(is.na(grass_025), NA, grass_025 + eps)  # preserve NA
ratio <- clamp(pasture_025 / denom, 0, 1)

drop_025 <- (ratio >= ALPHA)
mask_025 <- ifel(drop_025, 1L, 0L)  # 1=drop, 0=keep; NA propagated

## --- 5) Write 0.25° mask + 0.05° replica -------------------------------------
out025 <- file.path(out_dir, fname)
writeRaster(mask_025, out025, overwrite = TRUE,
            gdal = safe_gdal("LOG1S"), NAflag = 255)

mask_005 <- disagg(mask_025, fact = 5, method = "near")
if (!isTRUE(compareGeom(mask_005, ref005, stopOnError = FALSE))) {
  mask_005 <- resample(mask_005, ref005, method = "near")
}
out005 <- file.path(out_dir, sub("_0p25", "_0p05_rep", fname))
writeRaster(mask_005, out005, overwrite = TRUE,
            gdal = safe_gdal("LOG1S"), NAflag = 255)

## --- 6) Quicklooks (optional) -------------------------------------------------
outpng <- file.path(ql_dir, sub("\\.tif$", ".png", basename(out025)))
if (isTRUE(REMAKE_QL) || !file.exists(outpng)) {
  png(outpng, width = 1400, height = 700, res = 120)
  op <- par(mfrow = c(1, 3), mar = c(3, 3, 3, 6))
  plot(grass_025,   main = "Grass fraction (0.25°)",   col = pal_green,   zlim = c(0, 1))
  plot(pasture_025, main = "Pasture fraction (0.25°)", col = pal_green, zlim = c(0, 1))
  plot(mask_025,
       main = sprintf("Mask (P≥%s)", tok(ALPHA)),
       col = c("#f0f0f0", "#d73027"),
       breaks = c(-0.5, 0.5, 1.5),
       legend = FALSE, axes = TRUE, box = TRUE)
  legend("bottomleft", fill = c("#f0f0f0", "#d73027"),
         legend = c("0 keep", "1 drop"), bty = "n")
  par(op); dev.off()
}

gc()
cat(glue("
Wrote:
  - {out025}
  - {out005}
Rule: drop if pasture/grass≥{ALPHA}
Source grass={GRASS_SOURCE}; window={Y0}-{Y1}; semantics: 1=drop, 0=keep, 255=NA
"))

