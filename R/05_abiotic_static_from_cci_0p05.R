## =============================================================================
# 05_abiotic_static_from_cci.R — Build multi-year abiotic mask (0.05°)
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

ref005  <- rast(cfg$grids$grid_005$ref_raster)
cci_dir <- cfg$paths$cci_dir

TAU_WATER <- as.numeric(Sys.getenv("TAU_WATER", "0.05"))
TAU_ICE   <- as.numeric(Sys.getenv("TAU_ICE", "0.05"))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_abiotic")
ql_dir  <- file.path(out_dir, "quicklooks")

dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)

ESACCI <- cfg$esa_cci$classes

.as_int <- function(x) as.integer(unlist(x, use.names = FALSE))

vals_water  <- .as_int(ESACCI$water)
vals_ice    <- .as_int(ESACCI$snow_ice)
nodata_vals <- unique(c(.as_int(ESACCI$nodata), 255L))

files <- list.files(cci_dir, "*.tif$", full.names = TRUE)
years <- as.integer(str_extract(basename(files), "(19|20)\\d{2}"))

plan <- data.frame(file = files, year = years, stringsAsFactors = FALSE)
plan <- plan[!is.na(plan$year), ]

all_pW <- all_pI <- list()

for (i in seq_len(nrow(plan))) {

  out_tif <- file.path(
    out_dir,
    sprintf(
      "mask_static_abiotic_CCI_tauW%s_tauI%s_0p05.tif",
      gsub("\\.", "p", sprintf("%.2f", TAU_WATER)),
      gsub("\\.", "p", sprintf("%.2f", TAU_ICE))
    )
  )

  SKIP_EXISTING <- TRUE
  if (SKIP_EXISTING && file.exists(out_tif)) {
    message("✓ Abiotic mask already exists — skipping: ", out_tif)
    quit(save = "no")
  }

  r <- rast(plan$file[i])

  if (is.na(crs(r)))
    crs(r) <- crs(ref005)

  r[r %in% nodata_vals] <- NA

  bin <- function(rr, codes)
    classify(rr, cbind(codes, 1), others = 0)

  pW <- resample(bin(r, vals_water), ref005, method = "average")
  pI <- resample(bin(r, vals_ice),   ref005, method = "average")

  all_pW[[i]] <- pW
  all_pI[[i]] <- pI
}

mean_pW <- mean(rast(all_pW))
mean_pI <- mean(rast(all_pI))

abiotic <- (mean_pW >= TAU_WATER) | (mean_pI >= TAU_ICE)
abi_mask <- ifel(abiotic, 1L, 0L)
names(abi_mask) <- "abiotic_drop"

writeRaster(
  abi_mask,
  out_tif,
  overwrite = TRUE,
  NAflag = 255,
  datatype = "INT1U",
  gdal = c(
    "TILED=YES",
    "COMPRESS=DEFLATE",
    "PREDICTOR=2",
    "ZLEVEL=6",
    "BIGTIFF=IF_SAFER"
  )
)

# -------------------------------------------------------------------------
# Global + AOI quicklooks
# -------------------------------------------------------------------------
quicklook_mask_all_aois(
  mask    = abi_mask,
  title   = sprintf("Abiotic mask — All years (pW≥%.2f ∨ pI≥%.2f)", TAU_WATER, TAU_ICE),
  tag     = "abiotic_all_years",
  cfg     = cfg,
  ql_root = ql_dir,
  down    = 4L,
  include_global  = TRUE,
  drop_global_key = FALSE
)

# -------------------------------------------------------------------------
# AOI individual quicklooks
# -------------------------------------------------------------------------
aoi_out <- file.path(ql_dir, "aois")
dir.create(aoi_out, recursive = TRUE, showWarnings = FALSE)

aois <- cfg$aois
pal  <- hcl.colors(2, "batlow")

for (nm in names(aois)) {

  aoi <- aois[[nm]]
  ext_aoi <- ext(aoi$lon_min, aoi$lon_max, aoi$lat_min, aoi$lat_max)

  r_aoi <- crop(abi_mask, ext_aoi, snap = "near")

  ofile <- file.path(aoi_out, sprintf("abiotic_mask_%s.png", nm))

  png(ofile, width = 900, height = 900, res = 120)
  plot(r_aoi, col = pal, axes = TRUE, main = sprintf("Abiotic mask\nAOI: %s", nm))
  dev.off()
}

# -------------------------------------------------------------------------
# Final reporting
# -------------------------------------------------------------------------
p_drop <- global(ifel(abi_mask, 1, 0), "mean", na.rm = TRUE)[1, 1]

gc()

cat(
  glue(
    "\nWritten: {out_tif}\n",
    "Years: {min(plan$year)}–{max(plan$year)}\n",
    "Drop thresholds — water={TAU_WATER}, ice={TAU_ICE}\n",
    "Global drop fraction: {sprintf('%.4f', p_drop)}\n",
    "AOI quicklooks → {aoi_out}\n"
  )
)
