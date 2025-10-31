## =============================================================================
# 07_luh_masks.R — Build LUH pasture and rangeland share maps (0.25° + 0.05°)
#
# Purpose
#   Compute mean LUH pasture and rangeland area fractions over a given
#   multi-year window, outputting both native 0.25° fields and nearest-neighbor
#   0.05° replicas for masking or visualization purposes.
#
# Inputs
#   - LUH2 states.nc file (cfg$luh2$states_nc)
#     subdatasets: pasture, rangeland (or range)
#   - Reference grids:
#       * 0.25° grid: cfg$grids$grid_025$ref_raster
#       * 0.05° grid: cfg$grids$grid_005$ref_raster
#
# Outputs
#   - m_LUH_pasture_{Y0–Y1}_0p25.tif
#   - m_LUH_rangeland_{Y0–Y1}_0p25.tif
#   - m_LUH_pasture_{Y0–Y1}_0p05_rep.tif (optional)
#   - m_LUH_rangeland_{Y0–Y1}_0p05_rep.tif (optional)
#   - Quicklook PNGs for 0.25° and 0.05° replicas
#
# Environment variables
#   LUH_AVG_START, LUH_AVG_END (integer; averaging window)
#   WRITE_005_REPLICA (logical, default TRUE)
#
# Dependencies
#   Packages: terra, yaml, glue
#   Sourced helpers: utils.R, io.R, geom.R, viz.R, options.R
#
# Processing overview
#   1) Read LUH pasture and rangeland subdatasets from states.nc.
#   2) Subset to [Y0–Y1] and compute temporal means (0–1 fractions).
#   3) Align both layers to the reference 0.25° grid.
#   4) Optionally replicate to 0.05° grid by nearest disaggregation.
#   5) Write outputs and quicklooks; report year window and semantics.
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

Y0 <- as.integer(Sys.getenv("LUH_AVG_START", cfg$project$years$cci_start))
Y1 <- as.integer(Sys.getenv("LUH_AVG_END", cfg$project$years$cci_end))
WRITE_005 <- as.logical(Sys.getenv("WRITE_005_REPLICA", "TRUE"))

ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)

out_dir <- cfg$paths$masks_luh_dir
ql_dir  <- cfg$paths$masks_luh_quick_dir
dir.create(out_dir, TRUE, showWarnings = FALSE)
pal_green <- hcl.colors(64, "Greens", rev = TRUE)

# --- load LUH (0.25°) ----------------------------------------------------------
luh_nc <- cfg$luh2$states_nc
v_pas  <- cfg$luh2$variables$pasture
v_rng  <- cfg$luh2$variables$rangeland %||% cfg$luh2$variables$range
if (!file.exists(luh_nc))
  stop("LUH not found: ", luh_nc)

pas <- rast(luh_nc, subds = v_pas)
rng <- rast(luh_nc, subds = v_rng)

# subset years and average in time
ty  <- suppressWarnings(as.integer(time(pas)))
idx <- which(ty >= Y0 & ty <= Y1)
if (!length(idx))
  stop("No LUH steps in ", Y0, "-", Y1)

pas_mu <- mean(pas[[idx]], na.rm = TRUE)
rng_mu <- mean(rng[[idx]], na.rm = TRUE)

# shares at 0.25°
m025_pas <- clamp(pas_mu, 0, 1)
names(m025_pas) <- "pasture_share"
m025_rng <- clamp(rng_mu, 0, 1)
names(m025_rng) <- "rangeland_share"

# align to canonical 0.25° grid
if (!compareGeom(m025_pas, ref025, stopOnError = FALSE)) {
  m025_pas <- resample(m025_pas, ref025, method = "bilinear")
}
if (!compareGeom(m025_rng, ref025, stopOnError = FALSE)) {
  m025_rng <- resample(m025_rng, ref025, method = "bilinear")
}

# --- write 0.25° ---------------------------------------------------------------
tag <- glue("{Y0}-{Y1}")
f_pas_025 <- file.path(out_dir, glue("m_LUH_pasture_{tag}_0p25.tif"))
f_rng_025 <- file.path(out_dir, glue("m_LUH_rangeland_{tag}_0p25.tif"))

writeRaster(
  m025_pas,
  f_pas_025,
  overwrite = TRUE,
  gdal = gdal_wopt("FLT4S")$gdal,
  NAflag = -9999
)
writeRaster(
  m025_rng,
  f_rng_025,
  overwrite = TRUE,
  gdal = gdal_wopt("FLT4S")$gdal,
  NAflag = -9999
)

# --- 0.05° replicas (nearest replication) -------------------------------------
if (WRITE_005) {
  parent025    <- rast(cfg$grids$grid_025$ref_raster)
  parent_to_005 <- disagg(parent025, fact = 5, method = "near")  # index-safe 5×5 tiles

  m005_pas <- resample(m025_pas, parent_to_005, method = "near")
  m005_rng <- resample(m025_rng, parent_to_005, method = "near")

  # snap to canonical 0.05°
  if (!compareGeom(m005_pas, ref005, stopOnError = FALSE))
    m005_pas <- resample(m005_pas, ref005, method = "near")
  if (!compareGeom(m005_rng, ref005, stopOnError = FALSE))
    m005_rng <- resample(m005_rng, ref005, method = "near")

  f_pas_005 <- file.path(out_dir, glue("m_LUH_pasture_{tag}_0p05_rep.tif"))
  f_rng_005 <- file.path(out_dir, glue("m_LUH_rangeland_{tag}_0p05_rep.tif"))

  writeRaster(
    m005_pas,
    f_pas_005,
    overwrite = TRUE,
    gdal = gdal_wopt("FLT4S")$gdal,
    NAflag = -9999
  )
  writeRaster(
    m005_rng,
    f_rng_005,
    overwrite = TRUE,
    gdal = gdal_wopt("FLT4S")$gdal,
    NAflag = -9999
  )
}

# --- quicklooks ----------------------------------------------------------------
# (NB: ql_dir is overridden here to <out_dir>/quicklooks)
ql_dir <- file.path(out_dir, "quicklooks")
dir.create(ql_dir, TRUE, showWarnings = FALSE)
plot_ql <- function(r,
                    title,
                    file,
                    zlim = c(0, 1),
                    pal = pal_green) {
  png(file, 1400, 700, res = 120)
  op <- par(mar = c(3, 3, 3, 6))
  terra::plot(
    r,
    main = title,
    col = pal,
    zlim = zlim,
    axes = TRUE,
    legend = TRUE
  )
  mtext("0..1 fraction", 3, line = 0.7, cex = 0.9)
  par(op)
  dev.off()
}

plot_ql(m025_pas,
        sprintf("LUH pasture %s (0.25°)", tag),
        file.path(ql_dir, sprintf("pasture_%s_0p25.png", tag)))
plot_ql(m025_rng,
        sprintf("LUH rangeland %s (0.25°)", tag),
        file.path(ql_dir, sprintf("rangeland_%s_0p25.png", tag)))

if (WRITE_005) {
  plot_ql(
    m005_pas,
    sprintf("LUH pasture %s (0.05° replica)", tag),
    file.path(ql_dir, sprintf("pasture_%s_0p05_rep.png", tag))
  )
  plot_ql(
    m005_rng,
    sprintf("LUH rangeland %s (0.05° replica)", tag),
    file.path(ql_dir, sprintf("rangeland_%s_0p05_rep.png", tag))
  )
}
gc()
cat(
  glue(
    "
LUH masks written:
  0.25°:
    - {basename(f_pas_025)}
    - {basename(f_rng_025)}
  0.05° replicas: {if (WRITE_005) 'written' else 'skipped'}
Window: {Y0}-{Y1}
Semantics: shares in 0..1; pasture = LUH pasture; rangeland = LUH rangeland
"
  )
)
