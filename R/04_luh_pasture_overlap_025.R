## =============================================================================
# 04_luh_pasture_overlap_025.R — Build LUH pasture-overlap mask (0.25° + 0.05°)
#
# Purpose
#   Identify regions where LUH pasture fraction overlaps substantially with
#   satellite-derived grass cover. Produces a “drop” mask (1=drop, 0=keep)
#   indicating likely managed grassland based on grass–pasture ratios.
#
# Inputs
#   - LUH2 states.nc (cfg$luh2$states_nc; subds = "pastr")
#   - Grass fraction (0.05°) from GRASS_SOURCE = CCI | GLC_FRAC | GLC_TEMP
#   - Reference grids: cfg$grids$grid_005$ref_raster, cfg$grids$grid_025$ref_raster
#   - Cell area raster (0.05°): cfg$grids$grid_005$area_raster
#
# Outputs
#   - mask_luh_overlap_{SRC}_Gmin{G}_Pmin{P}_alpha{A}_{Y0–Y1}_0p25.tif
#   - mask_luh_overlap_{SRC}_Gmin{G}_Pmin{P}_alpha{A}_{Y0–Y1}_0p05_rep.tif
#   - Optional quicklook PNG
#
# Environment variables
#   GRASS_SOURCE (string, default "GLC_FRAC")
#   G_MIN        (numeric, default 0.05)
#   P_MIN        (numeric, default 0.20)
#   ALPHA        (numeric, default 0.50)
#   LUH_AVG_START, LUH_AVG_END (integer; year window)
#   REMAKE_QL    (logical, default TRUE)
#
# Dependencies
#   Packages: terra, yaml, glue, stringr
#   Sourced helpers: utils.R, io.R, geom.R, viz.R, options.R
#
# Processing overview
#   1) Derive mean grass fraction (0.05°) from chosen source over [Y0–Y1].
#   2) Aggregate to 0.25° using area-weighted mean.
#   3) Compute mean LUH pasture fraction (0.25°) for [Y0–Y1].
#   4) Drop where grass≥G_MIN, pasture≥P_MIN, and (pasture/grass)≥ALPHA.
#   5) Write 0.25° mask and 0.05° replica; optionally draw quicklook.
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(yaml)
  library(glue)
  library(stringr)
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

ref005  <- rast(cfg$grids$grid_005$ref_raster)
ref025  <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)

# --- params (env) --------------------------------------------------------------
GRASS_SOURCE <- toupper(Sys.getenv("GRASS_SOURCE", "GLC_FRAC"))  # CCI|GLC_FRAC
REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "TRUE"))
G_MIN  <- as.numeric(Sys.getenv("G_MIN", "0.05"))        # min grass to act
P_MIN  <- as.numeric(Sys.getenv("P_MIN", "0.20"))        # min pasture to act
ALPHA  <- as.numeric(Sys.getenv("ALPHA", "0.50"))        # pasture/grass threshold
Y0 <- get_int_env("LUH_AVG_START", cfg$project$years$cci_start)
Y1 <- get_int_env("LUH_AVG_END", cfg$project$years$cci_end)

stopifnot(is.finite(G_MIN),
          is.finite(P_MIN),
          is.finite(ALPHA),
          is.finite(Y0),
          is.finite(Y1))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_luh_overlap")
ql_dir  <- file.path(out_dir, "quicklooks")
dir.create(out_dir, TRUE, FALSE)
dir.create(ql_dir, TRUE, FALSE)



# --- 1) grass fraction @0.05° over [Y0..Y1] -----------------------------------
grass_005 <- switch(
  GRASS_SOURCE,
  "CCI" = {
    frac_dir <- cfg$paths$cci_out_dir
    f <- list.files(frac_dir, "ESACCI_frac_\\d{4}_0p05\\.tif$", full.names =
                      TRUE)
    if (!length(f))
      stop("No CCI fraction files in: ", frac_dir)
    yrs <- as.integer(str_extract(basename(f), "\\d{4}"))
    keep <- which(yrs >= Y0 & yrs <= Y1)
    if (!length(keep))
      stop("No CCI fraction years in window ", Y0, "-", Y1)
    stk <- rast(lapply(f[keep], function(x)
      rast(x)[["frac_grass"]]))
    mean(stk, na.rm = TRUE)
  },
  "GLC_FRAC" = {
    p <- file.path(cfg$paths$glc_out_dir, "GLC_frac_0p05.tif")
    if (!file.exists(p))
      stop("GLC_FRAC path not found: ", p)
    rast(p)
  },
  "GLC_TEMP" = {
    p <- file.path(cfg$paths$glc_out_dir, "glc_cat_yearstack_0p05.tif")
    if (!file.exists(p))
      stop("GLC yearstack not found: ", p)
    s <- rast(p)
    y <- suppressWarnings(as.integer(substr(names(s), 2, 5)))
    keep <- which(y >= Y0 & y <= Y1)
    if (!length(keep))
      stop("No GLC years in window ", Y0, "-", Y1)
    vec_int <- function(x)
      as.integer(unlist(x, use.names = FALSE))
    grass_vals <- vec_int(cfg$glc$classes$grassland)
    is_grass <- classify(s[[keep]], cbind(grass_vals, 1), others = 0)
    app(is_grass, mean, na.rm = TRUE)
  },
  stop(
    "Unknown GRASS_SOURCE='",
    GRASS_SOURCE,
    "'. Use CCI | GLC_FRAC | GLC_TEMP."
  )
)
if (!isTRUE(compareGeom(grass_005, ref005, stopOnError = FALSE)))
  grass_005 <- resample(grass_005, ref005, method = "bilinear")
names(grass_005) <- "grass"

# --- 2) aggregate grass 0.05° → 0.25° (area-weighted mean) --------------------
num <- aggregate(
  grass_005 * area005,
  fact = 5,
  fun = function(x)
    sum(x, na.rm = TRUE)
)
den <- aggregate((!is.na(grass_005)) * area005,
                 fact = 5,
                 fun = function(x)
                   sum(x, na.rm = TRUE)
)
grass_025 <- ifel(den == 0, NA, num / den)
if (!isTRUE(compareGeom(grass_025, ref025, stopOnError = FALSE)))
  grass_025 <- resample(grass_025, ref025, method = "near")
names(grass_025) <- "grass_025"

# --- 3) LUH pasture @0.25° (time mean) ---------------------------------------
luh_nc <- cfg$luh2$states_nc
if (!file.exists(luh_nc))
  stop("LUH file not found: ", luh_nc)
v_pas <- cfg$luh2$variables$pasture
# sanity: does subdataset exist?
sds_names <- try(names(sds(luh_nc)), silent = TRUE)
if (inherits(sds_names, "try-error"))
  sds_names <- character(0)
if (!(v_pas %in% sds_names))
  stop("LUH subdataset '",
       v_pas,
       "' not found. Available: ",
       paste(sds_names, collapse = ", "))

pas <- rast(luh_nc, subds = v_pas)
ty  <- suppressWarnings(as.integer(time(pas)))
if (all(is.na(ty)))
  stop("LUH time axis missing/NA for '", v_pas, "'. Check file.")
keep <- which(ty >= Y0 & ty <= Y1)
if (!length(keep))
  stop("No LUH timesteps in window ", Y0, "-", Y1)
pasture_025 <- mean(pas[[keep]], na.rm = TRUE)
if (!isTRUE(compareGeom(pasture_025, ref025, stopOnError = FALSE)))
  pasture_025 <- resample(pasture_025, ref025, method = "bilinear")
pasture_025 <- clamp(pasture_025, 0, 1)
names(pasture_025) <- "pasture"

# --- 4) decision at 0.25° -----------------------------------------------------
# Avoid base::pmax on rasters; add small epsilon safely
eps <- 1e-9
denom <- ifel(is.na(grass_025), NA, grass_025 + eps)  # keep NA where grass is NA
ratio <- clamp(pasture_025 / denom, 0, 1)
drop_025 <- (grass_025 >= G_MIN) &
  (pasture_025 >= P_MIN) & (ratio >= ALPHA)
mask_025 <- ifel(drop_025, 1L, 0L)  # 1=drop, 0=keep, NA propagated

# --- 5) write 0.25° + 0.05° replica ------------------------------------------

out025 <- file.path(out_dir, fname)
writeRaster(
  mask_025,
  out025,
  overwrite = TRUE,
  gdal = gdal_wopt("LOG1S")$gdal,
  NAflag = 255
)

mask_005 <- disagg(mask_025, fact = 5, method = "near")
if (!isTRUE(compareGeom(mask_005, ref005, stopOnError = FALSE)))
  mask_005 <- resample(mask_005, ref005, method = "near")
out005 <- file.path(out_dir, sub("_0p25", "_0p05_rep", fname))
writeRaster(
  mask_005,
  out005,
  overwrite = TRUE,
  gdal = gdal_wopt("LOG1S")$gdal,
  NAflag = 255
)

# --- 6) quicklooks -------------------------------------------------------------
if (REMAKE_QL || !file.exists(...)) {
  png(file.path(ql_dir, sub("\\.tif$", ".png", basename(out025))), 1400, 700, res =
        120)
  op <- par(mfrow = c(1, 3), mar = c(3, 3, 3, 6))

  terra::plot(grass_025,
              main = "Grass fraction (0.25°)",
              col = pal_grass,
              zlim = c(0, 1))
  terra::plot(pasture_025,
              main = "Pasture fraction (0.25°)",
              col = pal_pasture,
              zlim = c(0, 1))
  terra::plot(
    mask_025,
    main = sprintf("Mask (α=%s; G≥%s; P≥%s)", tok(ALPHA), tok(G_MIN), tok(P_MIN)),
    col = c("#f0f0f0", "#d73027"),
    breaks = c(-0.5, 0.5, 1.5),
    legend = FALSE,
    axes = TRUE,
    box = TRUE
  )
  legend(
    "bottomleft",
    fill = c("#f0f0f0", "#d73027"),
    legend = c("0 keep", "1 drop"),
    bty = "n"
  )
  par(op)
  dev.off()
}
gc()
cat(
  glue(
    "
Wrote:
  - {out025}
  - {out005}
Rule: drop if grass≥{G_MIN} & pasture≥{P_MIN} & pasture/grass≥{ALPHA}
Source grass={GRASS_SOURCE}; window={Y0}-{Y1}; semantics: 1=drop, 0=keep, 255=NA
"
  )
)
