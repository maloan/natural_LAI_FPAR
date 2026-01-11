## =============================================================================
# 05_glc_stack_0p05.R — Build annual GLC_FCS30D categorical yearstack (0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "viz.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.9)

ref005 <- rast(cfg$grids$grid_005$ref_raster)

glc_dir <- cfg$paths$glc_dir
out_dir <- cfg$paths$glc_out_dir
ql_dir <- file.path(out_dir, "quicklooks")
stack_out <- file.path(out_dir, "glc_cat_yearstack_0p05.tif")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
OVERWRITE <- as.logical(Sys.getenv("OVERWRITE", "FALSE"))
REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))

years_wanted <- as.integer(cfg$glc$years)
classes <- cfg$glc$classes

vec_int <- function(x) as.integer(unlist(x, use.names = FALSE))
cropland_vals <- vec_int(classes$cropland)
urban_vals <- vec_int(classes$urban)
nodata_vals <- vec_int(classes$nodata)

years_for_ql <- intersect(c(1990, 2000, 2010, 2020), years_wanted)

glc_quicklook_layers <- function(cat_r) {
  rC <- classify(cat_r, cbind(cropland_vals, 1), others = 0)
  rU <- classify(cat_r, cbind(urban_vals, 1), others = 0)
  rast(list(frac_cropland = rC, frac_urban = rU))
}

# --- discover files (base R) --------------------------------------------------
all_files <- list.files(glc_dir, pattern = "\\.tif$", full.names = TRUE)
if (!length(all_files)) stop("No GLC GeoTIFFs found in: ", glc_dir)

get_year <- function(p) {
  m <- regexpr("(19|20)\\d{2}", basename(p), perl = TRUE)
  if (m[1] < 0) {
    return(NA_integer_)
  }
  as.integer(substr(basename(p), m[1], m[1] + attr(m, "match.length") - 1))
}

years <- vapply(all_files, get_year, integer(1))
ok <- !is.na(years) & years %in% years_wanted
paths <- all_files[ok]
years <- years[ok]

ord <- order(years)
paths <- paths[ord]
years <- years[ord]

stopifnot(length(paths) > 0)

message(
  "Found ", length(paths),
  " GLC rasters; span [", min(years), "..", max(years), "]."
)

# --- if stack exists: optionally quicklooks and exit --------------------------
if (file.exists(stack_out) && SKIP_EXISTING && !OVERWRITE) {
  message("✓ Yearstack exists — skipping rebuild: ", stack_out)

  if (REMAKE_QL && length(years_for_ql)) {
    s <- rast(stack_out)
    for (yr in years_for_ql) {
      nm <- sprintf("Y%04d", yr)
      if (nm %in% names(s)) {
        ql_layers <- glc_quicklook_layers(s[[nm]])
        quicklook_all_aois(
          frac = ql_layers,
          year = yr,
          cfg = cfg,
          ql_root = ql_dir,
          down = 4L,
          include_global = TRUE,
          drop_global_key = TRUE
        )
      }
    }
  }
  quit(save = "no")
}

# --- rebuild stack ------------------------------------------------------------
bands <- vector("list", length(paths))

for (i in seq_along(paths)) {
  yr <- years[i]
  f <- paths[i]

  message("→ Processing ", basename(f))

  r <- rast(f)[[1]]
  if (is.na(crs(r))) crs(r) <- crs(ref005)

  if (!compareGeom(r, ref005, stopOnError = FALSE)) {
    r <- resample(r, ref005, method = "near")
  }

  if (length(nodata_vals)) {
    r <- terra::subst(r, nodata_vals, NA)
  }

  names(r) <- sprintf("Y%04d", yr)
  bands[[i]] <- r

  if (yr %in% years_for_ql) {
    ql_layers <- glc_quicklook_layers(r)
    quicklook_all_aois(
      frac = ql_layers,
      year = yr,
      cfg = cfg,
      ql_root = ql_dir,
      down = 4L,
      include_global = TRUE,
      drop_global_key = TRUE
    )
  }

  gc()
}
# wopt <- wopt_int(opts$SPEED_OVER_SIZE)
stack <- rast(bands)
writeRaster(stack, stack_out,
  overwrite = TRUE
)
gc()

cat("Wrote GLC yearstack: ", stack_out, "\n", sep = "")
