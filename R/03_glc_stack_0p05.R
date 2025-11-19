## =============================================================================
# 03_glc_stack_0p05.R — Build annual GLC_FCS30D categorical yearstack (0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(yaml)
  library(dplyr)
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

ref005 <- rast(cfg$grids$grid_005$ref_raster)

# --- paths ---------------------------------------------------------------------
glc_dir   <- cfg$paths$glc_dir
out_dir   <- cfg$paths$glc_out_dir
ql_dir    <- file.path(out_dir, "quicklooks")
stack_out <- file.path(out_dir, "glc_cat_yearstack_0p05.tif")

dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)

# toggles
SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
OVERWRITE     <- as.logical(Sys.getenv("OVERWRITE", "FALSE"))
REMAKE_QL     <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))

# --- config & class lists -----------------------------------------------------
vec_int <- function(x) as.integer(unlist(x, use.names = FALSE))

years_wanted <- as.integer(cfg$glc$years)
classes <- cfg$glc$classes

cropland_vals <- vec_int(classes$cropland)
urban_vals    <- vec_int(classes$urban)
nodata_vals   <- vec_int(classes$nodata)

# ----------------- Quicklook helper -------------------------------------------
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

  rC <- classify(cat_r, cbind(cropland_vals, 1), others = 0)
  rU <- classify(cat_r, cbind(urban_vals, 1), others = 0)
  mask <- rast(list(mask_cropland = rC, mask_urban = rU))
  names(mask) <- c("mask_cropland", "mask_urban")
  return(mask)
}

# --- discover per-year files ---------------------------------------------------
files <- list.files(glc_dir, pattern = "\\.tif$", full.names = TRUE) |>
  tibble(path = _) |>
  mutate(year = as.integer(str_extract(basename(path), "(19|20)\\d{2}"))) |>
  filter(!is.na(year), year %in% years_wanted) |>
  arrange(year)

stopifnot(nrow(files) > 0)

missing <- setdiff(years_wanted, files$year)
if (length(missing)) {
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

HAVE_STACK <- file.exists(stack_out)
years_for_ql <- intersect(c(1990, 2000, 2010, 2020), years_wanted)

if (SKIP_EXISTING && HAVE_STACK && !OVERWRITE) {

  message("✓ Yearstack exists — skipping rebuild: ", stack_out)

  if (REMAKE_QL && length(years_for_ql)) {
    s <- rast(stack_out)
    for (yr in years_for_ql) {
      nm <- sprintf("Y%04d", yr)
      if (!nm %in% names(s)) next
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

  quit(save = "no")
}

# --- rebuild stack --------------------------------------------------------------
bands <- list()

for (i in seq_len(nrow(files))) {

  yr <- files$year[i]
  f  <- files$path[i]
  r  <- rast(f)

  if (is.na(crs(r))) crs(r) <- crs(ref005)

  seen <- unique(values(r))
  seen <- seen[is.finite(seen)]

  allowed <- unique(unlist(classes[c(
    "cropland","urban","grassland","bare","water","snow_ice"
  )], use.names = FALSE))

  extra <- setdiff(seen, allowed)
  if (length(extra)) {
    warning("GLC unexpected codes in ",
            basename(f),
            ": ",
            paste(head(extra, 20), collapse = ", "))
  }

  if (!compareGeom(r, ref005, stopOnError = FALSE)) {
    message("Resampling ", basename(f), " → 0.05° grid (near)")
    r <- resample(r, ref005, method = "near")
  }

  if (length(nodata_vals))
    r[r %in% nodata_vals] <- NA

  names(r) <- sprintf("Y%04d", yr)
  bands[[length(bands) + 1]] <- r

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

stack <- rast(bands)

writeRaster(
  stack,
  stack_out,
  overwrite = TRUE,
  wopt = wopt_int(opts$SPEED_OVER_SIZE)
)

cat("Wrote glc yearstack: ", stack_out, "\n", sep = "")
gc()
