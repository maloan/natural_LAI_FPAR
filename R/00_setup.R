# =============================================================================
# 00_setup.R — Initialize config, reference grids, areas (0.05°/0.25°)
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(terra)
  library(here)
})

source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))
terraOptions(progress = 1, memfrac = 0.25)

epsg4326 <- "EPSG:4326"
run_tag <- Sys.getenv("run_tag", unset = "alpha_0.1")
alpha_cci <- as.numeric(Sys.getenv("alpha_cci", "0.1"))
GDAL_OPTS <- gdal_co_int()

#  paths
in_dirs <- list(
  lai_nc_dir = exp_(here("data-raw", "LAI", "lai_1982-2024")),
  fpar_nc_dir = exp_(here("data-raw", "FPAR", "fpar_1982-2024")),
  cci_dir = exp_(here("data-raw", "ESACCI", "ESACCI_1992-2020")),
  luh2_dir = exp_(here("data-raw", "LUH2_v2h")),
  glc_dir = exp_(here("data-raw", "GLC_FCS30D")),
  valid_tiles_info = exp_(here(
    "src", "valid_tiles_info_0p05_full_10deg.rds"
  )),
  bilinear_ref = exp_(here("src", "refgrid_0p05.nc"))
)

dirs <- list(
  ref_dir = here("src"),
  out_root = here("output", run_tag),
  eval_dir = here("output", run_tag, "eval"),
  georef_lai_0p05_dir = here("data", "georef", "georef_lai_0p05"),
  georef_fpar_0p05_dir = here("data", "georef", "georef_fpar_0p05"),
  cci_out_dir = here("data", "frac", "cci_frac_0p05"),
  cci_quick_dir = here("data", "frac", "cci_frac_0p05", "quicklooks"),
  glc_out_dir = here("data", "frac", "glc_frac_0p05"),
  glc_quick_dir = here("data", "frac", "glc_frac_0p05", "quicklooks"),
  masks_root_dir = here("output", run_tag, "masks"),
  masks_cci_dir = here("output", run_tag, "masks", "mask_cci"),
  masks_cci_quick_dir = here("output", run_tag, "masks", "mask_cci", "quicklooks"),
  masks_glc_dir = here("output", run_tag, "masks", "mask_glc"),
  masks_glc_quick_dir = here("output", run_tag, "masks", "mask_glc", "quicklooks"),
  masks_luh_dir = here("output", run_tag, "masks", "mask_luh"),
  masks_luh_quick_dir = here("output", run_tag, "masks", "mask_luh", "quicklooks"),
  mask_nonvegetated_dir = here("output", run_tag, "masks", "mask_nonvegetated")
)

# masked dirs
for (v in c("LAI", "FPAR")) {
  for (m in c("CCI", "GLC")) {
    for (r in c("0p05", "0p25", "0p5")) {
      key <- sprintf("masked_%s_%s_%s_dir", tolower(v), tolower(m), r)
      dirs[[key]] <- here(
        "output",
        run_tag,
        paste0("masked_", r),
        v,
        sprintf("masked_%s_%s", v, m)
      )
    }
  }
}

invisible(lapply(
  unique(unlist(dirs)),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

#  product paths
resolutions <- list(
  `0p05` = list(nrc = c(3600, 7200), deg = 0.05),
  `0p25` = list(nrc = c(720, 1440), deg = 0.25),
  `0p5`  = list(nrc = c(360, 720), deg = 0.5)
)

ref <- lapply(names(resolutions), function(res) {
  # Build reference rasters and netcdf files for the specified resolution.
  # Returns a list with paths to the tif and nc files, as well as the nrc and
  # deg values.
  list(
    tif = here("src", sprintf("ref_%s.tif", res)),
    nc = here("src", sprintf("ref_%s.nc", res)),
    nrc = resolutions[[res]]$nrc,
    deg = resolutions[[res]]$deg
  )
})
names(ref) <- names(resolutions)

area <- lapply(names(resolutions), function(res) {
  # Build area rasters and netcdf files for the specified resolution. Returns a
  # list with paths to the tif and nc files, as well as the deg value.
  list(
    tif = here("src", sprintf("area_%s_km2.tif", res)),
    nc = here("src", sprintf("area_%s_km2.nc", res)),
    deg = resolutions[[res]]$deg
  )
})
names(area) <- names(resolutions)

#  config (write once per run_tag)
cfg_path <- here("config", sprintf("config_%s.yml", run_tag))
dir.create(dirname(cfg_path),
  recursive = TRUE,
  showWarnings = FALSE
)

cfg <- list(
  project = list(
    run_tag = run_tag,
    alpha_cci = alpha_cci,
    name = "SNU_LAI_FPAR_natmask_global",
    crs = epsg4326,
    years = list(
      lai_start = 1982,
      lai_end = 2024,
      cci_start = 1992,
      cci_end = 2020,
      glc_start = 1992,
      glc_end = 2020
    )
  )
)

cfg$paths <- utils::modifyList(in_dirs, dirs)

cfg$grids <- list(
  grid_005 = list(
    ref_raster = ref$`0p05`$nc,
    area_raster = area$`0p05`$nc
  ),
  grid_025 = list(
    ref_raster = ref$`0p25`$nc,
    area_raster = area$`0p25`$nc
  ),
  grid_05 = list(
    ref_raster = ref$`0p5`$nc,
    area_raster = area$`0p5`$nc
  )
)

cfg$variables <- list(
  produce = list("LAI", "FPAR"),
  lai = list(
    nc_var_name_primary = "LAI",
    nc_var_name_fallback = "auto_first_variable",
    units = "m2 m-2",
    nc_lon_name = "lon",
    nc_lat_name = "lat",
    nc_time_names = c("time", "time_counter")
  ),
  fpar = list(
    nc_var_name_primary = "FPAR",
    nc_var_name_fallback = "auto_first_variable",
    units = "1",
    nc_lon_name = "lon",
    nc_lat_name = "lat",
    nc_time_names = c("time", "time_counter")
  )
)

cfg$resampling <- list(categorical = "near", continuous = "bilinear")

cfg$esa_cci <- list(
  version = "v2.0.7",
  classes = list(
    nodata = 0,
    cropland = c(10, 11, 12, 20),
    cls30 = 30,
    cls40 = 40,
    grassland = 130,
    urban = 190,
    water = 210,
    snow_ice = 220
  ),
  mask_window_years = c(1992, 2020),
  clean_majority_threshold = 0.5,
  clean_operator = "<=",
  weights = list(cls30 = 0.75, cls40 = 0.25)
)

cfg$glc <- list(
  product = "GLC_FCS30D",
  tiles_dir = cfg$paths$glc_dir,
  nodata_in = 0,
  classes = list(
    cropland = c(10, 11, 12, 20),
    grassland = 130,
    urban = 190,
    bare = c(200, 201, 202),
    water = 210,
    snow_ice = 220,
    nodata = c(0, 250)
  ),
  years = c(1985, 1990, 1995, 2000:2022),
  mask_window_years = c(1992, 2020),
  clean_majority_threshold = 0.5,
  clean_operator = "<="
)

cfg$luh2 <- list(
  states_nc = here("data-raw", "LUH2_v2h", "states.nc"),
  variables = list(
    cropland_components = c("c3ann", "c4ann", "c3per", "c4per", "c3nfx"),
    pasture = "pastr",
    pasture_range = "range",
    urban = "urban"
  )
)

cfg$thresholds <- list(
  cu_fraction_max_025 = c(0.03, 0.05, 0.10, 0.20),
  baseline_T = 0.05
)

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

#  build reference rasters + netcdf
for (k in names(ref)) {
  tif <- ref[[k]]$tif
  nc <- ref[[k]]$nc
  nrc <- ref[[k]]$nrc
  deg <- ref[[k]]$deg


  r <- rast(
    nrows = nrc[1],
    ncols = nrc[2],
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 90,
    crs = epsg4326
  )
  values(r) <- NA_real_
  writeRaster(r, tif, overwrite = TRUE, gdal = GDAL_OPTS)
  terra::writeCDF(r, nc, overwrite = TRUE)
}

#  build area rasters + netcdf
for (k in names(area)) {
  tif <- area[[k]]$tif
  nc <- area[[k]]$nc
  deg <- area[[k]]$deg

  a <- cellSize(rast(ref[[k]]$tif), unit = "km")
  writeRaster(a, tif, overwrite = TRUE, gdal = GDAL_OPTS)
  terra::writeCDF(a, nc, overwrite = TRUE)
}
# sanity checks & manifest
ref_rasts <- lapply(ref, function(r) {
  rast(r$tif)
})
area_rasts <- lapply(area, function(a) {
  rast(a$tif)
})

man <- data.frame(
  file = sapply(ref, function(r) {
    basename(r$tif)
  }),
  rows = sapply(ref_rasts, nrow),
  cols = sapply(ref_rasts, ncol),
  res_x = sapply(ref_rasts, function(r) {
    res(r)[1]
  }),
  res_y = sapply(ref_rasts, function(r) {
    res(r)[2]
  }),
  xmin = sapply(ref_rasts, xmin),
  xmax = sapply(ref_rasts, xmax),
  ymin = sapply(ref_rasts, ymin),
  ymax = sapply(ref_rasts, ymax),
  crs = sapply(ref_rasts, function(r) {
    crs(r)
  }),
  global_area_km2 = sapply(area_rasts, function(a) {
    as.numeric(global(a, "sum", na.rm = TRUE)[1, 1])
  }),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

ea <- man$global_area_km2
stopifnot(abs(ea[1] - ea[2]) < 1e-2 * ea[1])
stopifnot(abs(mean(ea) - 510072000) < 2e6)

write.csv(man, here("src", "manifest_00.csv"), row.names = FALSE)

# Clean up temporary objects
rm(ref_rasts, area_rasts, resolutions, GDAL_OPTS)
gc()
