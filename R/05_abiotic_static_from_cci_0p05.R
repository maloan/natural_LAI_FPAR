# 05_abiotic_static_from_cci_all_years.R
# Build a multi-year abiotic mask using ESA-CCI categorical rasters (~300 m)

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
source(file.path(ROOT, "R", "00_utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))
source(file.path(ROOT, "R", "viz.R"))
cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
cci_dir <- cfg$paths$cci_dir

TAU_WATER <- as.numeric(Sys.getenv("TAU_WATER", "0.05"))
TAU_ICE   <- as.numeric(Sys.getenv("TAU_ICE",   "0.05"))
TAU_BARE  <- as.numeric(Sys.getenv("TAU_BARE",  "0.30"))

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_abiotic")
ql_dir  <- file.path(out_dir, "quicklooks")
dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)

ESACCI <- cfg$esa_cci$classes
.as_int <- function(x) as.integer(unlist(x, use.names = FALSE))
vals_water <- .as_int(ESACCI$water)
vals_ice   <- .as_int(ESACCI$snow_ice)
vals_bare  <- .as_int(ESACCI$bare)
nodata_vals <- unique(c(.as_int(ESACCI$nodata), 255L))

files <- list.files(cci_dir, "\.tif$", full.names = TRUE)
years <- as.integer(str_extract(basename(files), "(19|20)\\d{2}"))
plan <- data.frame(file = files, year = years, stringsAsFactors = FALSE)
plan <- plan[!is.na(plan$year), ]

all_pW <- all_pI <- all_pB <- list()

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
if (SKIP_EXISTING && file.exists(out_tif)) {
  message("✓ Abiotic mask already exists — skipping: ", out_tif)
  quit(save = "no")
}


for (i in seq_len(nrow(plan))) {
  r <- rast(plan$file[i])
  if (is.na(crs(r))) crs(r) <- crs(ref005)
  r[r %in% nodata_vals] <- NA

  bin <- function(rr, codes) classify(rr, cbind(codes, 1), others = 0)
  # Convert binary presence/absence to fraction of pixel-year presence
  pW <- resample(bin(r, vals_water), ref005, method = "average")
  pI <- resample(bin(r, vals_ice),   ref005, method = "average")
  pB <- resample(bin(r, vals_bare),  ref005, method = "average")

  all_pW[[i]] <- pW
  all_pI[[i]] <- pI
  all_pB[[i]] <- pB
}

mean_pW <- mean(rast(all_pW))
mean_pI <- mean(rast(all_pI))
mean_pB <- mean(rast(all_pB))

abiotic <- (mean_pW >= TAU_WATER) | (mean_pI >= TAU_ICE) | (mean_pB >= TAU_BARE)
abi_mask <- ifel(abiotic, 1L, 0L)
names(abi_mask) <- "abiotic_drop"

out_tif <- file.path(out_dir, sprintf("mask_static_abiotic_CCI_all_tauW%s_tauI%s_tauB%s_0p05.tif",
                                      gsub("\\.", "p", sprintf("%.2f", TAU_WATER)),
                                      gsub("\\.", "p", sprintf("%.2f", TAU_ICE)),
                                      gsub("\\.", "p", sprintf("%.2f", TAU_BARE))))

writeRaster(abi_mask, out_tif, overwrite = TRUE, NAflag = 255,
            datatype = "INT1U",
            gdal = c("TILED=YES","COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=6","BIGTIFF=IF_SAFER"))

quicklook_mask_all_aois(
  mask    = abi_mask,
  title   = sprintf("Abiotic mask — All years (pW≥%.2f ∨ pI≥%.2f ∨ pB≥%.2f)", TAU_WATER, TAU_ICE, TAU_BARE),
  tag     = "abiotic_all_years",
  cfg     = cfg,
  ql_root = ql_dir,
  down    = 4L,
  include_global  = TRUE,
  drop_global_key = FALSE
)

p_drop <- global(ifel(abi_mask, 1, 0), "mean", na.rm = TRUE)[1,1]
cat(glue("\nWritten: {out_tif}\nYears: {min(plan$year)}–{max(plan$year)}\n",
         "Drop thresholds — water={TAU_WATER}, ice={TAU_ICE}, bare={TAU_BARE}\n",
         "Global drop fraction: {sprintf('%.4f', p_drop)}\n"))
