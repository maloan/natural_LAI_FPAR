## =============================================================================
# 09_luh_masks.R — Build LUH pasture and rangeland share maps (0.25° + 0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "viz.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.25)

Y0 <- as.integer(Sys.getenv("LUH_AVG_START", cfg$project$years$cci_start))
Y1 <- as.integer(Sys.getenv("LUH_AVG_END", cfg$project$years$cci_end))
WRITE_005 <- as.logical(Sys.getenv("WRITE_005", "TRUE"))

ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)

out_dir <- cfg$paths$masks_luh_dir
ql_dir <- file.path(out_dir, "quicklooks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

luh_nc <- cfg$luh2$states_nc
v_pas <- cfg$luh2$variables$pasture
v_rng <- cfg$luh2$variables$rangeland %||% cfg$luh2$variables$range
if (!file.exists(luh_nc)) stop("LUH not found: ", luh_nc, call. = FALSE)

pas <- rast(luh_nc, subds = v_pas)
rng <- rast(luh_nc, subds = v_rng)
wopt <- gdal_wopt("FLT4S")
ty <- suppressWarnings(as.integer(time(pas)))
idx <- which(ty >= Y0 & ty <= Y1)
if (!length(idx)) stop("No LUH steps in ", Y0, "-", Y1, call. = FALSE)

m025_pas <- clamp(mean(pas[[idx]], na.rm = TRUE), 0, 1)
m025_rng <- clamp(mean(rng[[idx]], na.rm = TRUE), 0, 1)
names(m025_pas) <- "pasture_share"
names(m025_rng) <- "rangeland_share"

# align once to canonical 0.25° grid
if (!compareGeom(m025_pas, ref025, stopOnError = FALSE)) {
  m025_pas <- resample(m025_pas, ref025, method = "bilinear")
}
if (!compareGeom(m025_rng, ref025, stopOnError = FALSE)) {
  m025_rng <- resample(m025_rng, ref025, method = "bilinear")
}

tag <- sprintf("%d-%d", Y0, Y1)

f_pas_025 <- file.path(out_dir, sprintf("m_LUH_pasture_%s_0p25.tif", tag))
f_rng_025 <- file.path(out_dir, sprintf("m_LUH_rangeland_%s_0p25.tif", tag))

writeRaster(m025_pas, f_pas_025,
  overwrite = TRUE,
  gdal = gdal_wopt("FLT4S"), NAflag = -9999
)
writeRaster(m025_rng, f_rng_025,
  overwrite = TRUE,
  gdal = gdal_wopt("FLT4S"), NAflag = -9999
)

# 0.05° replicas (nearest-neighbour replication)
if (WRITE_005) {
  m005_pas <- disagg(m025_pas, fact = 5, method = "near")
  m005_rng <- disagg(m025_rng, fact = 5, method = "near")

  if (!compareGeom(m005_pas, ref005, stopOnError = FALSE)) {
    m005_pas <- resample(m005_pas, ref005, method = "near")
  }
  if (!compareGeom(m005_rng, ref005, stopOnError = FALSE)) {
    m005_rng <- resample(m005_rng, ref005, method = "near")
  }

  f_pas_005 <- file.path(
    out_dir,
    sprintf("m_LUH_pasture_%s_0p05_rep.tif", tag)
  )
  f_rng_005 <- file.path(
    out_dir,
    sprintf("m_LUH_rangeland_%s_0p05_rep.tif", tag)
  )

  writeRaster(m005_pas, f_pas_005,
    overwrite = TRUE,
    gdal = gdal_wopt("FLT4S"), NAflag = -9999
  )
  writeRaster(m005_rng, f_rng_005,
    overwrite = TRUE,
    gdal = gdal_wopt("FLT4S"), NAflag = -9999
  )
}

# quicklooks (optional but cheap)
write_quicklook_raster(m025_pas, file.path(
  ql_dir,
  sprintf("pasture_%s_0p25.png", tag)
),
title = sprintf("LUH pasture %s (0.25°)", tag)
)
write_quicklook_raster(m025_rng, file.path(
  ql_dir, sprintf("rangeland_%s_0p25.png", tag)
),
title = sprintf("LUH rangeland %s (0.25°)", tag)
)

if (WRITE_005) {
  write_quicklook_raster(m005_pas, file.path(
    ql_dir, sprintf("pasture_%s_0p05_rep.png", tag)
  ),
  title = sprintf("LUH pasture %s (0.05° replica)", tag)
  )
  write_quicklook_raster(m005_rng, file.path(
    ql_dir, sprintf("rangeland_%s_0p05_rep.png", tag)
  ),
  title = sprintf("LUH rangeland %s (0.05° replica)", tag)
  )
}

gc()
cat(sprintf(
  "LUH masks written (window %s):\n  0.25°: %s, %s\n  0.05° replicas: %s\n",
  tag, basename(f_pas_025), basename(f_rng_025), if (WRITE_005) {
    "written"
  } else {
    "skipped"
  }
))
