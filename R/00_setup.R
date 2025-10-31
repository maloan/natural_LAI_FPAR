## =============================================================================
# 00_setup.R — Initialize config, reference grids, areas, and AOIs (0.05° / 0.25°)
#
# Purpose
#   Bootstrap the SNU_LAI_FPAR environment by generating reference grids,
#   global cell-area rasters, AOI masks, and a structured config.yml that
#   defines paths, datasets, and processing parameters for the workflow.
#
# Inputs
#   - Environment variable SNU_LAI_FPAR_ROOT (default: ~/GitHub/natural_LAI_FPAR)
#   - Raw data directories: LAI, FPAR, ESA-CCI, LUH2, GLC_FCS30D
#   - Utility scripts: R/utils.R, R/io.R
#
# Outputs
#   - config/config.yml — unified configuration and paths
#   - src/ref_0p05.tif, ref_0p25.tif — reference grids
#   - src/area_0p05_km2.tif, area_0p25_km2.tif — cell areas (km²)
#   - src/aoi_0p05.tif, aoi_0p25.tif — AOI masks
#   - src/manifest_00.csv — metadata and validation summary
#
# Environment variables
#   RECOMPUTE_REFS  (logical, default FALSE) — rebuild reference grids
#   RECOMPUTE_AREAS (logical, default FALSE) — rebuild area rasters
#   RECOMPUTE_AOIS  (logical, default FALSE) — rebuild AOI masks
#   SILENT_TIMING   (logical, default FALSE)
#   GDAL_TILED      (logical, default TRUE)
#   RUN_TAG         (string, default "baseline")
#
# Processing overview
#   1) Define directory structure and dataset paths.
#   2) Write project config (grids, variables, thresholds, masks).
#   3) Generate and validate reference rasters (0.05°, 0.25°).
#   4) Compute cell areas and AOI masks.
#   5) Sanity-check CRS, extents, and total Earth area.
#   6) Write manifest for reproducibility.
## =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(terra)
})

# --- helpers & root ------------------------------------------------------------
ROOT <- normalizePath(path.expand(
  Sys.getenv("SNU_LAI_FPAR_ROOT", "~/GitHub/natural_LAI_FPAR")
),
winslash = "/",
mustWork = FALSE)
source(file.path(ROOT, "R", "utils.R"))
source(file.path(ROOT, "R", "io.R"))

terraOptions(progress = 1, memfrac = 0.25)

# --- toggles -------------------------------------------------------------------
RECOMPUTE_REFS  <- as.logical(Sys.getenv("RECOMPUTE_REFS", "FALSE"))
RECOMPUTE_AREAS <- as.logical(Sys.getenv("RECOMPUTE_AREAS", "FALSE"))
RECOMPUTE_AOIS  <- as.logical(Sys.getenv("RECOMPUTE_AOIS", "FALSE"))
SILENT_TIMING   <- as.logical(Sys.getenv("SILENT_TIMING", "FALSE"))
GDAL_TILED      <- as.logical(Sys.getenv("GDAL_TILED", "TRUE"))

EPSG4326 <- "EPSG:4326"

# Centralized GTiff creation options (Float32 by default), with tiling toggle
GDAL_OPTS <- {
  co <- gdal_wopt("FLT4S")$gdal
  if (!GDAL_TILED)
    co <- setdiff(co, "TILED=YES")
  co
}

# --- dirs ----------------------------------------------------------------------
REF_DIR <- file.path(ROOT, "src")
RUN_TAG <- Sys.getenv("RUN_TAG", unset = "baseline")
OUT_DIR <- file.path(ROOT, "output", RUN_TAG)
COMMON_OUT_DIR <- file.path(ROOT, "data")

# snapshot of inputs (expandables kept explicit so config.yml is self-explanatory)
in_dirs <- list(
  lai_nc_dir       = exp_(file.path(ROOT, "data-raw/LAI/lai_1982-2021")),
  fpar_nc_dir      = exp_(file.path(ROOT, "data-raw/FPAR/fpar_1982-2021")),
  cci_dir          = exp_(file.path(
    ROOT, "data-raw/ESACCI/ESACCI_1992-2020"
  )),
  luh2_dir         = exp_(file.path(ROOT, "data-raw/LUH2_v2h")),
  glc_dir          = exp_(file.path(ROOT, "data-raw/GLC_FCS30D")),
  valid_tiles_info = exp_(file.path(
    ROOT, "src/valid_tiles_info_0p05_full_10deg.rds"
  )),
  bilinear_ref     = exp_(file.path(ROOT, "src/refgrid_0p05.nc"))
)

# canonical output tree
dirs <- list(
  ref_dir   = REF_DIR,
  out_root  = OUT_DIR,

  eval_dir = file.path(OUT_DIR, "eval"),

  # georef (0.05°)
  georef_lai_0p05_dir  = file.path(COMMON_OUT_DIR, "georef", "georef_lai_0p05"),
  georef_fpar_0p05_dir = file.path(COMMON_OUT_DIR, "georef", "georef_fpar_0p05"),

  # fractions + quicklooks
  cci_out_dir     = file.path(COMMON_OUT_DIR, "frac", "cci_frac_0p05"),
  cci_quick_dir   = file.path(COMMON_OUT_DIR, "frac", "cci_frac_0p05", "quicklooks"),
  glc_out_dir     = file.path(COMMON_OUT_DIR, "frac", "GLC_frac_0p05"),
  glc_quick_dir   = file.path(COMMON_OUT_DIR, "frac", "GLC_frac_0p05", "quicklooks"),

  # masks + quicklooks
  masks_root_dir      = file.path(OUT_DIR, "masks"),
  masks_cci_dir       = file.path(OUT_DIR, "masks", "mask_cci"),
  masks_cci_quick_dir = file.path(OUT_DIR, "masks", "mask_cci", "quicklooks"),
  masks_glc_dir       = file.path(OUT_DIR, "masks", "mask_glc"),
  masks_glc_quick_dir = file.path(OUT_DIR, "masks", "mask_glc", "quicklooks"),
  masks_luh_dir       = file.path(OUT_DIR, "masks", "mask_luh"),
  masks_luh_quick_dir = file.path(OUT_DIR, "masks", "mask_luh", "quicklooks"),

  # masked 0.05°
  masked_lai_cci_005_dir  = file.path(
    OUT_DIR,
    "masked_0p05",
    "LAI",
    "masked_LAI_CCI",
    "masked_lai_0p05"
  ),
  masked_lai_glc_005_dir  = file.path(
    OUT_DIR,
    "masked_0p05",
    "LAI",
    "masked_LAI_GLC",
    "masked_lai_0p05"
  ),
  masked_fpar_cci_005_dir = file.path(
    OUT_DIR,
    "masked_0p05",
    "FPAR",
    "masked_FPAR_CCI",
    "masked_fpar_0p05"
  ),
  masked_fpar_glc_005_dir = file.path(
    OUT_DIR,
    "masked_0p05",
    "FPAR",
    "masked_FPAR_GLC",
    "masked_fpar_0p05"
  ),

  # masked 0.25°
  masked_lai_cci_025_dir  = file.path(
    OUT_DIR,
    "masked_0p25",
    "LAI",
    "masked_LAI_CCI",
    "masked_lai_0p25"
  ),
  masked_lai_glc_025_dir  = file.path(
    OUT_DIR,
    "masked_0p25",
    "LAI",
    "masked_LAI_GLC",
    "masked_lai_0p25"
  ),
  masked_fpar_cci_025_dir = file.path(
    OUT_DIR,
    "masked_0p25",
    "FPAR",
    "masked_FPAR_CCI",
    "masked_fpar_0p25"
  ),
  masked_fpar_glc_025_dir = file.path(
    OUT_DIR,
    "masked_0p25",
    "FPAR",
    "masked_FPAR_GLC",
    "masked_fpar_0p25"
  )
)

# ensure tree exists
invisible(lapply(
  unique(unlist(dirs)),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# --- write config.yml ----------------------------------------------------------
cfg_path <- file.path(ROOT, "config", "config.yml")


cfg <- list()

cfg$project <- list(
  run_tag = RUN_TAG,
  name = "SNU_LAI_FPAR_natmask_global",
  crs  = EPSG4326,
  aoi_extent = list(
    lon_min = 20,
    lon_max = 140,
    lat_min = 35,
    lat_max = 70
  ),
  years = list(
    lai_start = 1992,
    lai_end = 2020,
    cci_start = 1992,
    cci_end = 2020,
    glc_start = 1985,
    glc_end = 2020
  )
)

# AOIs for quicklooks (named extents)
cfg$aois <- cfg$project$aois <- list(
  india       = list(
    lon_min =  68,
    lon_max =  98,
    lat_min =  6,
    lat_max = 37
  ),
  china       = list(
    lon_min =  73,
    lon_max = 135,
    lat_min = 18,
    lat_max = 54
  ),
  south_korea = list(
    lon_min = 124,
    lon_max = 131,
    lat_min = 33,
    lat_max = 40
  ),
  switzerland = list(
    lon_min =   6,
    lon_max =  11,
    lat_min = 46,
    lat_max = 48
  ),
  usa         = list(
    lon_min = -170,
    lon_max = -65,
    lat_min = 18,
    lat_max = 72
  ),
  europe      = list(
    lon_min = -25,
    lon_max =  45,
    lat_min = 33,
    lat_max = 72
  ),
  se_asia     = list(
    lon_min =  92,
    lon_max = 120,
    lat_min = -12,
    lat_max = 23
  ),
  amazon      = list(
    lon_min = -80,
    lon_max = -45,
    lat_min = -20,
    lat_max =  6
  ),
  australia_se = list(
    lon_min = 134,
    lon_max = 154,
    lat_min = -40,
    lat_max = -24
  ),
  indonesia   = list(
    lon_min =  95,
    lon_max = 141,
    lat_min = -11,
    lat_max =  6
  ),
  nile_delta  = list(
    lon_min =  28,
    lon_max =  33,
    lat_min = 29,
    lat_max = 32
  ),
  pampas      = list(
    lon_min = -66,
    lon_max = -57,
    lat_min = -39,
    lat_max = -30
  ),
  corn_belt   = list(
    lon_min = -105,
    lon_max = -80,
    lat_min = 35,
    lat_max = 49
  ),
  uk          = list(
    lon_min = -10,
    lon_max =   3,
    lat_min = 49,
    lat_max = 60
  ),
  west_africa = list(
    lon_min = -20,
    lon_max =  20,
    lat_min =  0,
    lat_max = 20
  ),
  global      = list(
    lon_min = -180,
    lon_max = 180,
    lat_min = -90,
    lat_max = 90
  )
)

cfg$paths <- utils::modifyList(in_dirs, dirs)

cfg$grids <- list(
  grid_005 = list(
    ref_raster = file.path(REF_DIR, "ref_0p05.tif"),
    resolution_deg = 0.05,
    extent = c(-180, 180, -90, 90)
  ),
  grid_025 = list(
    ref_raster = file.path(REF_DIR, "ref_0p25.tif"),
    resolution_deg = 0.25,
    extent = c(-180, 180, -90, 90)
  )
)

cfg$variables <- list(
  produce = list("LAI", "FPAR"),
  lai = list(
    nc_var_name_primary  = "LAI",
    nc_var_name_fallback = "auto_first_variable",
    units = "m2 m-2",
    clamp = list(min = 0, max = 8),
    nc_lon_name   = "lon",
    nc_lat_name   = "lat",
    nc_time_names = c("time", "time_counter")
  ),
  fpar = list(
    nc_var_name_primary  = "FPAR",
    nc_var_name_fallback = "auto_first_variable",
    units = "1",
    clamp = list(min = 0, max = 1),
    nc_lon_name   = "lon",
    nc_lat_name   = "lat",
    nc_time_names = c("time", "time_counter")
  )
)

cfg$resampling <- list(categorical = "near", continuous = "bilinear")

cfg$esa_cci <- list(
  version = "v2.0.7",
  classes = list(
    nodata   = 0,
    cropland = c(10, 11, 12, 20),
    cls30    = 30,
    cls40    = 40,
    grassland = 130,
    urban    = 190,
    bare     = c(200, 201, 202),
    water    = 210,
    snow_ice = 220
  ),
  mask_window_years = c(1992, 2020),
  clean_majority_threshold = 0.5,
  clean_operator = "<=",
  weights = list(cls30 = 0.75, cls40 = 0.25)
)

cfg$glc <- list(
  product   = "GLC_FCS30D",
  tiles_dir = file.path(ROOT, "data-raw/GLC_FCS30D"),
  nodata_in = 0,
  classes   = list(
    cropland = c(10, 11, 12, 20),
    grassland = 130,
    urban    = 190,
    bare     = c(200, 201, 202),
    water    = 210,
    snow_ice = 220,
    nodata   = c(0, 250)
  ),
  years = c(1985, 1990, 1995, 2000:2022),
  mask_window_years = c(1985, 2020),
  clean_majority_threshold = 0.5,
  clean_operator = "<="
)

cfg$luh2 <- list(
  states_nc = file.path(in_dirs$luh2_dir, "states.nc"),
  variables = list(
    cropland_components = c("c3ann", "c4ann", "c3per", "c4per", "c3nfx"),
    pasture             = "pastr",
    pasture_range       = "range",
    urban               = "urban"
  )
)

cfg$thresholds <- list(cu_fraction_max_025 = c(0.03, 0.05, 0.10, 0.20),
                       baseline_T = 0.05)

cfg$naming <- list(
  masked005      = "LAI_{pipeline}_{mask_key}_{yyyymm}_0p05_masked.tif",
  masked025      = "LAI_{pipeline}_{mask_key}_{yyyymm}_0p25_masked.tif",
  final025       = "LAI_{pipeline}_{mask_key}_{yyyymm}_0p25_luh2_T={T}.tif",
  T_token_format = "0p%02d",
  masked005_fpar = "FPAR_{pipeline}_{mask_key}_{yyyymm}_0p05_masked.tif",
  masked025_fpar = "FPAR_{pipeline}_{mask_key}_{yyyymm}_0p25_masked.tif",
  final025_fpar  = "FPAR_{pipeline}_{mask_key}_{yyyymm}_0p25_luh2_T={T}.tif"
)

yaml::write_yaml(cfg, cfg_path)

# --- reference rasters ---------------------------------------------------------
ref005_path <- file.path(ROOT, "src", "ref_0p05.tif")
ref025_path <- file.path(ROOT, "src", "ref_0p25.tif")

make_ref <- function(nrows, ncols, crs) {
  r <- rast(
    nrows = nrows,
    ncols = ncols,
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 90,
    crs = crs
  )
  values(r) <- NA_real_
  r
}

if (RECOMPUTE_REFS || !file.exists(ref005_path)) {
  writeRaster(
    make_ref(3600, 7200, EPSG4326),
    ref005_path,
    overwrite = TRUE,
    gdal = GDAL_OPTS,
    description = gtiff_desc(
      product = "ref_grid",
      grid_deg = 0.05,
      extras = list(empty = TRUE)
    )
  )
}
if (RECOMPUTE_REFS || !file.exists(ref025_path)) {
  writeRaster(
    make_ref(720, 1440, EPSG4326),
    ref025_path,
    overwrite = TRUE,
    gdal = GDAL_OPTS,
    description = gtiff_desc(
      product = "ref_grid",
      grid_deg = 0.25,
      extras = list(empty = TRUE)
    )
  )
}

# --- area rasters --------------------------------------------------------------
area005_path <- file.path(REF_DIR, "area_0p05_km2.tif")
area025_path <- file.path(REF_DIR, "area_0p25_km2.tif")
aoi005_path  <- file.path(REF_DIR, "aoi_0p05.tif")
aoi025_path  <- file.path(REF_DIR, "aoi_0p25.tif")


if (RECOMPUTE_AREAS || !file.exists(area005_path)) {
  a005 <- cellSize(rast(ref005_path), unit = "km")
  writeRaster(
    a005,
    area005_path,
    overwrite = TRUE,
    gdal = GDAL_OPTS,
    description = gtiff_desc(
      product = "cell_area_km2",
      grid_deg = 0.05,
      extras = list(source = "terra::cellSize(unit='km')")
    )
  )
}
if (RECOMPUTE_AREAS || !file.exists(area025_path)) {
  a025 <- cellSize(rast(ref025_path), unit = "km")
  writeRaster(
    a025,
    area025_path,
    overwrite = TRUE,
    gdal = GDAL_OPTS,
    description = gtiff_desc(
      product = "cell_area_km2",
      grid_deg = 0.25,
      extras = list(source = "terra::cellSize(unit='km')")
    )
  )
}

# --- AOI masks -----------------------------------------------------------------
AOI <- cfg$project$aoi_extent
poly <- as.polygons(ext(AOI$lon_min, AOI$lon_max, AOI$lat_min, AOI$lat_max), crs = EPSG4326)
aoi_str <- sprintf("lon[%g..%g],lat[%g..%g]",
                   AOI$lon_min,
                   AOI$lon_max,
                   AOI$lat_min,
                   AOI$lat_max)

if (RECOMPUTE_AOIS || !file.exists(aoi005_path)) {
  a005 <- rasterize(poly, rast(ref005_path), field = 1)
  writeRaster(
    a005,
    aoi005_path,
    overwrite = TRUE,
    gdal = GDAL_OPTS,
    description = gtiff_desc(
      product = "aoi_mask",
      grid_deg = 0.05,
      extras = list(
        aoi = aoi_str,
        inside = 1,
        outside = "NA"
      )
    )
  )
}
if (RECOMPUTE_AOIS || !file.exists(aoi025_path)) {
  a025 <- rasterize(poly, rast(ref025_path), field = 1)
  writeRaster(
    a025,
    aoi025_path,
    overwrite = TRUE,
    gdal = GDAL_OPTS,
    description = gtiff_desc(
      product = "aoi_mask",
      grid_deg = 0.25,
      extras = list(
        aoi = aoi_str,
        inside = 1,
        outside = "NA"
      )
    )
  )
}

# --- finalize paths in config --------------------------------------------------
cfg <- cfg_read()
cfg$grids$grid_005$area_raster <- area005_path
cfg$grids$grid_025$area_raster <- area025_path
cfg$grids$grid_005$aoi_raster  <- aoi005_path
cfg$grids$grid_025$aoi_raster  <- aoi025_path
yaml::write_yaml(cfg, file.path(ROOT, "config", "config.yml"))

# --- manifest & sanity checks --------------------------------------------------
ref005 <- rast(ref005_path)
ref025 <- rast(ref025_path)
area005 <- rast(area005_path)
area025 <- rast(area025_path)
aoi005  <- rast(aoi005_path)
aoi025  <- rast(aoi025_path)

is_int_like <- function(x, tol = 1e-10)
  abs(x - round(x)) < tol
is_snapped  <- function(r,
                        origin_x = -180,
                        origin_y = -90,
                        tol = 1e-10) {
  rx <- res(r)[1]
  ry <- res(r)[2]
  qx <- (xmin(r) - origin_x) / rx
  qy <- (ymin(r) - origin_y) / ry
  is_int_like(qx, tol) && is_int_like(qy, tol)
}

stopifnot(
  terra::same.crs(ref005, ref025),
  isTRUE(all.equal(res(ref005) * 5, res(ref025))),
  isTRUE(all.equal(ext(ref005), ext(-180, 180, -90, 90))),
  isTRUE(all.equal(ext(ref025), ext(-180, 180, -90, 90))),
  is_snapped(ref005),
  is_snapped(ref025)
)

man <- data.frame(
  file = c(basename(ref005_path), basename(ref025_path)),
  rows = c(nrow(ref005), nrow(ref025)),
  cols = c(ncol(ref005), ncol(ref025)),
  res_x = c(res(ref005)[1], res(ref025)[1]),
  res_y = c(res(ref005)[2], res(ref025)[2]),
  xmin = c(xmin(ref005), xmin(ref025)),
  xmax = c(xmax(ref005), xmax(ref025)),
  ymin = c(ymin(ref005), ymin(ref025)),
  ymax = c(ymax(ref005), ymax(ref025)),
  crs  = c(crs(ref005, describe = TRUE), crs(ref025, describe = TRUE)),
  global_area_km2 = c(as.numeric(global(
    area005, "sum", na.rm = TRUE
  )[1, 1]), as.numeric(global(
    area025, "sum", na.rm = TRUE
  )[1, 1])),
  aoi_area_km2 = c(as.numeric(global(
    mask(area005, aoi005), "sum", na.rm = TRUE
  )[1, 1]), as.numeric(global(
    mask(area025, aoi025), "sum", na.rm = TRUE
  )[1, 1])),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

ea <- man$global_area_km2
stopifnot(abs(ea[1] - ea[2]) < 1e-2 * ea[1])       # cross-grid area near-identical
stopifnot(abs(mean(ea) - 510072000) < 2e6)         # ~0.4% tolerance to Earth area

write.csv(man, file.path(REF_DIR, "manifest_00.csv"), row.names = FALSE)

gc()
cat("Done 00_setup.R\n")
