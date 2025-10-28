## =============================================================================
# 04_glc_used_0p05.R — Build “used ≥ N years” mask from GLC_FCS30D yearstack
#
# Purpose
#   Derive a binary “used-land” mask (1=drop, 0=keep) from the annual
#   GLC_FCS30D categorical yearstack by counting how often a pixel was classified
#   as cropland or urban and applying a persistence threshold (≥N years).
#
# Inputs
#   - Yearstack: glc_cat_yearstack_0p05.tif (from 03_glc_stack_0p05.R)
#   - Reference grid: cfg$grids$grid_005$ref_raster
#   - Class codes: cfg$glc$classes (cropland, urban, nodata)
#   - Year window: cfg$project$years$cci_start : cci_end
#
# Outputs
#   - mask_used_ge{N}_{Y1–Y2}_0p05.tif (Byte; 1=drop, 0=keep, NA=255)
#   - glc_counts_crop_urban_0p05.tif (Int; per-pixel counts)
#   - Quicklook PNGs for global and AOIs
#
# Environment variables
#   USED_N_YEARS  (integer, default 3)
#   SKIP_EXISTING (logical, default FALSE)
#   OVERWRITE     (logical, default TRUE)
#   REMAKE_QL     (logical, default TRUE)
#
# Dependencies
#   Packages: terra, yaml, glue
#   Sourced helpers: utils.R, io.R, geom.R, viz.R, options.R
#
# Processing overview
#   1) Read GLC yearstack and align to 0.05° template grid.
#   2) Count years where pixel ∈ {cropland, urban}.
#   3) Mark pixels used in ≥N years → drop mask.
#   4) Write mask and count rasters; generate global/AOI quicklooks.
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(yaml)
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

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "FALSE"))
OVERWRITE     <- as.logical(Sys.getenv("OVERWRITE", "TRUE"))
REMAKE_QL     <- as.logical(Sys.getenv("REMAKE_QL", "TRUE"))

# --- paths ---------------------------------------------------------------------
glc_out_dir <- path.expand(cfg$paths$glc_out_dir)
masks_dir <- path.expand(cfg$paths$masks_glc_dir)
ql_dir      <- file.path(masks_dir, "quicklooks")
dir.create(masks_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)

# --- config / class codes ------------------------------------------------------
N_YEARS <- as.integer(Sys.getenv("USED_N_YEARS", "3"))
if (!is.finite(N_YEARS) ||
    N_YEARS < 1)
  stop("USED_N_YEARS must be a positive integer")

vec_int      <- function(x)
  as.integer(unlist(x, use.names = FALSE))
classes      <- cfg$glc$classes
cropland_vals <- vec_int(classes$cropland)
urban_vals    <- vec_int(classes$urban)
nodata_vals   <- vec_int(classes$nodata)

# --- template & year window ----------------------------------------------------
tmpl        <- rast(cfg$grids$grid_005$ref_raster)
cci_win     <- cfg$project$years$cci_start:cfg$project$years$cci_end
Y1 <- min(cci_win)
Y2 <- max(cci_win)

# --- find yearstack ------------------------------------------------------------
stack <- file.path(glc_out_dir, "glc_cat_yearstack_0p05.tif")
stack_path <- stack[file.exists(stack)]
if (!length(stack_path))
  stop("GLC yearstack not found under: ", glc_out_dir)
stack_path <- stack_path[1]
cat("Using GLC yearstack: ", stack_path, "\n", sep = "")

# --- read & align --------------------------------------------------------------
s <- rast(stack_path)                 # VRT or GTiff
if (is.na(crs(s)))
  crs(s) <- crs(tmpl)
if (!compareGeom(s, tmpl, stopOnError = FALSE)) {
  message("Resampling yearstack → 0.05° template (near)")
  s <- resample(s, tmpl, method = "near")
}
# sanitize nodata
if (length(nodata_vals))
  s[s %in% nodata_vals] <- NA

# --- restrict to CCI window (once) --------------------------------------------
keep_idx <- which(substr(names(s), 2, 5) %in% as.character(cci_win))
s <- s[[keep_idx]]
cat("Year window for used-counts: ",
    Y1,
    "..",
    Y2,
    " (",
    nlyr(s),
    " years)\n",
    sep = "")

# --- cores ---------------------------------------------------------------------
ncores <- opts$N_WORKERS
if (is.null(ncores) ||
    is.na(ncores) || !is.finite(ncores))
  ncores <- 1L
ncores <- max(1L, as.integer(ncores))

# --- per-pixel counts & mask ---------------------------------------------------
cnt_cropland <- app(
  s,
  fun  = function(v, vals)
    sum(v %in% vals, na.rm = TRUE),
  vals = cropland_vals,
  cores = ncores
)
cnt_urban <- app(
  s,
  fun  = function(v, vals)
    sum(v %in% vals, na.rm = TRUE),
  vals = urban_vals,
  cores = ncores
)
names(cnt_cropland) <- "cnt_cropland"
names(cnt_urban)    <- "cnt_urban"

cnt_total <- cnt_cropland + cnt_urban
used_geN  <- (cnt_total >= N_YEARS)
used_byte <- ifel(used_geN, 1L, 0L)  # 1=drop, 0=keep

# --- outputs -------------------------------------------------------------------
out_used   <- file.path(masks_dir,
                        sprintf("mask_used_ge%d_%d-%d_0p05.tif", N_YEARS, Y1, Y2))
out_counts <- file.path(masks_dir, "glc_counts_crop_urban_0p05.tif")

if (!(SKIP_EXISTING && file.exists(out_used)) || OVERWRITE) {
  writeRaster(
    used_byte,
    out_used,
    overwrite = TRUE,
    wopt = wopt_byte(opts$SPEED_OVER_SIZE, na = 255L)
  )
}
if (!(SKIP_EXISTING && file.exists(out_counts)) || OVERWRITE) {
  writeRaster(
    c(cnt_cropland, cnt_urban),
    out_counts,
    overwrite = TRUE,
    wopt = wopt_int(opts$SPEED_OVER_SIZE)
  )
}

# --- quicklooks ---------------------------------------------------------------
tag <- sprintf("glc_used_ge%d_%d-%d", N_YEARS, Y1, Y2)
ql_probe <- file.path(ql_dir,
                      "global",
                      sprintf("quicklook_mask_global_%s.png", tag))

if (REMAKE_QL || !file.exists(ql_probe)) {
  quicklook_mask_all_aois(
    mask    = used_byte,
    title   = sprintf("GLC mask (used \u2265 %d years)", N_YEARS),
    tag     = tag,
    cfg     = cfg,
    ql_root = ql_dir,
    down    = 4L,
    include_global  = TRUE,
    drop_global_key = FALSE
  )
}

# --- summary ------------------------------------------------------------------
p_used_global <- tryCatch(
  global(ifel(used_byte, 1, 0), "mean", na.rm = TRUE)[1, 1],
  error = function(e)
    NA_real_
)

cat(
  glue(
    "
Wrote:
  - {out_used}
  - {out_counts}
Quicklooks in: {ql_dir}
Used≥{N_YEARS} proportion (global): {sprintf('%.4f', p_used_global)}
"
  )
)
gc()
cat("Done\n")
