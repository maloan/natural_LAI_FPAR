## =============================================================================
# 03_glc_stack_0p05.R — Build annual GLC_FCS30D categorical yearstack (0.05°)
#
# Purpose
#   Assemble annual categorical GLC_FCS30D land-cover rasters (1985–2022) into
#   a single multi-band GeoTIFF aligned to the 0.05° global reference grid.
#   Generates periodic quicklooks for cropland and urban masks.
#
# Inputs
#   - Annual GLC_FCS30D GeoTIFFs: cfg$paths$glc_dir
#   - Reference grid: cfg$grids$grid_005$ref_raster (EPSG:4326)
#   - Class codes and years: cfg$glc$classes, cfg$glc$years
#
# Outputs
#   - glc_cat_yearstack_0p05.tif — multi-layer categorical stack (1 band/year)
#   - Quicklooks: cropland/urban binary maps for selected years (e.g., 1990–2020)
#
# Environment variables
#   SKIP_EXISTING (logical, default TRUE)  — skip if stack exists
#   OVERWRITE     (logical, default FALSE) — force rebuild even if exists
#   REMAKE_QL     (logical, default FALSE) — regenerate quicklooks only
#
# Dependencies
#   Packages: terra, yaml, dplyr, stringr, glue
#   Sourced helpers: utils.R, io.R, geom.R, viz.R, options.R
#
# Processing overview
#   1) Identify all annual categorical rasters within config year range.
#   2) Validate class codes and enforce CRS/alignment to 0.05° grid.
#   3) Stack yearly layers into a single multi-band raster.
#   4) Write GeoTIFF and generate global/AOI quicklooks (cropland, urban).
## =============================================================================

# Load packages
suppressPackageStartupMessages({
  library(terra)
  library(yaml)
  library(dplyr)
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
source(file.path(ROOT, "R", "utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))
source(file.path(ROOT, "R", "viz.R"))
source(file.path(ROOT, "R", "options.R"))
cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)
ref005 <- rast(path.expand(cfg$grids$grid_005$ref_raster))

# --- paths ---------------------------------------------------------------------
glc_dir <- path.expand(cfg$paths$glc_dir)        # where annual categorical rasters live
out_dir     <- path.expand(cfg$paths$glc_out_dir)    # where stack will be written
ql_dir      <- file.path(out_dir, "quicklooks")
stack_out <- file.path(out_dir, "glc_cat_yearstack_0p05.tif")

dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)

# toggles
SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
OVERWRITE     <- as.logical(Sys.getenv("OVERWRITE", "FALSE"))
REMAKE_QL  <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))
# --- config & class lists -----------------------------------------------------
vec_int <- function(x) {
  # create vector of integers
  as.integer(unlist(x, use.names = FALSE))
}
# load years from config
years_wanted <- as.integer(cfg$glc$years)
classes <- cfg$glc$classes
# convert class ID lists from config to integer vectors
cropland_vals <- vec_int(classes$cropland)
urban_vals    <- vec_int(classes$urban)
nodata_vals   <- vec_int(classes$nodata)

# ----------------- Quicklook helper -------------------------------------------
# Build 0/1 display layers and name them like the fraction quicklook expects
glc_quicklook_layers <- function(cat_r, cropland_vals, urban_vals) {
  rC <- classify(cat_r, cbind(cropland_vals, 1), others = 0)
  rU <- classify(cat_r, cbind(urban_vals, 1), others = 0)
  x  <- rast(list(frac_cropland = rC, frac_urban = rU))
  names(x) <- c("frac_cropland", "frac_urban")
  x
}

# ----------------- Plotting helper --------------------------------------------
plot_quicklooks <- function(cat_r,
                            cropland_vals,
                            urban_vals,
                            yr,
                            cfg,
                            ql_dir) {
  # make a binary cropland mask (1 where class ∈ cropland_vals, else 0).
  rC <- classify(cat_r, cbind(cropland_vals, 1), others = 0)
  # make a binary urban mask (1 where class ∈ urban_vals, else 0).
  rU <- classify(cat_r, cbind(urban_vals, 1), others = 0)
  # stack the two masks into a 2-layer raster.
  mask <- rast(list(mask_cropland = rC, mask_urban = rU))
  # set layer names for quicklooks
  names(mask) <- c("mask_cropland", "mask_urban")
  return(mask)
}

# --- discover per-year files ---------------------------------------------------
# Expect filenames that contain the 4-digit year, e.g. "..._2000_0p05.tif"
files <- list.files(glc_dir, pattern = "\\.tif$", full.names = TRUE) |>
  tibble(path = _) |>
  # define year from base name
  mutate(year = as.integer(str_extract(basename(path), "(19|20)\\d{2}"))) |>
  filter(!is.na(year), year %in% years_wanted) |>
  arrange(year)

# stop if number of input files is null
stopifnot(nrow(files) > 0)
# calculate number of years missing in input
missing <- setdiff(years_wanted, files$year)
if (length(missing)) {
  # print warning if input years are missing
  warning("Missing glc years: ", paste(missing, collapse = ", "))
}

message(
  "Found ",
  nrow(files),
  " glc annual rasters; span [",
  min(files$year),
  "..",
  max(files$year),
  "]."
)

HAVE_STACK <- file.exists(stack_out) # Check if stack already exists
years_for_ql <- intersect(c(1990, 2000, 2010, 2020), years_wanted)

if (SKIP_EXISTING && HAVE_STACK && !OVERWRITE) {
  # If stack exists and no overwrite, do not rebuild
  message("✓ Yearstack exists — skipping rebuild: ", stack_out)
  if (REMAKE_QL && length(years_for_ql)) {
    s <- rast(stack_out)
    for (yr in years_for_ql) {
      nm <- sprintf("Y%04d", yr)
      if (!nm %in% names(s))
        next
      cat_r <- s[[nm]]
      ql_layers <- glc_quicklook_layers(cat_r, cropland_vals, urban_vals)
      quicklook_all_aois(
        frac    = ql_layers,
        year    = yr,
        cfg     = cfg,
        ql_root = ql_dir,
        down    = 4L,
        include_global  = TRUE,
        drop_global_key = TRUE
      )
    }
  }

  # when the existing stack is kept exit the script immediately
  quit(save = "no")
} else {
  # rebuild the yearstack
  bands <- list() # container for per-year rasters
  for (i in seq_len(nrow(files))) {
    # for each file in the input files
    yr <- files$year[i]
    f  <- files$path[i]
    r  <- rast(f) # open the year’s categorical raster
    if (is.na(crs(r))) {
      # if CRS missing, set it to match ref005
      crs(r) <- crs(ref005)
    }
    seen <- unique(values(r))
    seen <- seen[is.finite(seen)]
    allowed <- unique(unlist(classes[c("cropland",
                                       "urban",
                                       "grassland",
                                       "bare",
                                       "water",
                                       "snow_ice")], use.names = FALSE))
    extra <- setdiff(seen, allowed)
    if (length(extra)) {
      warning("GLC unexpected codes in ",
              basename(f),
              ": ",
              paste(head(extra, 20), collapse = ", "))
    }

    if (!compareGeom(r, ref005, stopOnError = FALSE)) {
      message("Resampling ", basename(f), " → 0.05° grid (near)")
      # if geometry doesn’t match ref005, resample with nearest neighbor
      r <- resample(r, ref005, method = "near")
    }
    if (length(nodata_vals)) {
      # Set nodata to NA
      r[r %in% nodata_vals] <- NA
    }
    # name the layer by year.
    names(r) <- sprintf("Y%04d", yr)
    # append to the list for later stacking
    bands[[length(bands) + 1]] <- r

    # a couple of quicklooks (optional): only for some years to keep it light
    # Same as above
    if (yr %in% c(1990, 2000, 2010, 2020)) {
      ql_layers <- glc_quicklook_layers(r, cropland_vals, urban_vals)
      quicklook_all_aois(
        frac    = ql_layers,
        year    = yr,
        cfg     = cfg,
        ql_root = ql_dir,
        down    = 4L,
        include_global  = TRUE,
        drop_global_key = TRUE
      )
    }
  }
  # combine the per-year rasters in bands into a multi-layer
  stack <- rast(bands)
  # save that stack to stack_out as a GeoTIFF,
  writeRaster(stack,
              stack_out,
              overwrite = TRUE,
              wopt = wopt_int(opts$SPEED_OVER_SIZE))
  cat("Wrote glc yearstack: ", stack_out, "\n", sep = "")
}
gc()
