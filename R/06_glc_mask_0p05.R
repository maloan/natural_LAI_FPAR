## =============================================================================
# 06_glc_mask_0p05.R — Build “used ≥ N years” mask from GLC_FCS30D yearstack
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "viz.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.9)

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
OVERWRITE <- as.logical(Sys.getenv("OVERWRITE", "FALSE"))

N_YEARS <- as.integer(Sys.getenv("USED_N_YEARS", "3"))

glc_out_dir <- cfg$paths$glc_out_dir
masks_dir <- cfg$paths$masks_glc_dir
dir.create(masks_dir, recursive = TRUE, showWarnings = FALSE)

tmpl <- rast(cfg$grids$grid_005$ref_raster)

classes <- cfg$glc$classes
cropland_vals <- as.integer(unlist(classes$cropland, use.names = FALSE))
urban_vals <- as.integer(unlist(classes$urban, use.names = FALSE))
nodata_vals <- as.integer(unlist(classes$nodata, use.names = FALSE))

# analysis window (match CCI window)
yrs <- cfg$project$years$cci_start:cfg$project$years$cci_end
Y1 <- min(yrs)
Y2 <- max(yrs)

stack_path <- file.path(glc_out_dir, "glc_cat_yearstack_0p05.tif")
if (!file.exists(stack_path)) {
  stop("GLC yearstack not found: ", stack_path)
}
message("Using GLC yearstack: ", stack_path)

s <- rast(stack_path)
if (is.na(crs(s))) {
  crs(s) <- crs(tmpl)
}
if (!compareGeom(s, tmpl, stopOnError = FALSE)) {
  s <- resample(s, tmpl, method = "near")
}
if (length(nodata_vals)) {
  s[s %in% nodata_vals] <- NA
}

year_from_name <- function(nm) as.integer(substr(nm, 2, 5)) # expects "Y1990"
layer_years <- vapply(names(s), year_from_name, integer(1))
keep_idx <- which(layer_years %in% yrs)
s <- s[[keep_idx]]
stopifnot(nlyr(s) > 0)

message(sprintf(
  "Year window for used-counts: %d..%d (%d years)",
  Y1, Y2, nlyr(s)
))

ncores <- max(1L, as.integer(opts$N_WORKERS %||% 1L))

cnt_cropland <- app(s, function(v, vals) sum(v %in% vals, na.rm = TRUE),
  vals = cropland_vals, cores = ncores
)
cnt_urban <- app(s, function(v, vals) sum(v %in% vals, na.rm = TRUE),
  vals = urban_vals, cores = ncores
)

names(cnt_cropland) <- "cnt_cropland"
names(cnt_urban) <- "cnt_urban"

cnt_total <- cnt_cropland + cnt_urban
used_byte <- ifel(cnt_total >= N_YEARS, 1L, 0L)

out_used <- file.path(
  masks_dir,
  sprintf(
    "mask_used_ge%d_%d-%d_0p05.tif",
    N_YEARS, Y1, Y2
  )
)
out_counts <- file.path(masks_dir, "glc_counts_crop_urban_0p05.tif")

if (OVERWRITE || !(SKIP_EXISTING && file.exists(out_used))) {
  writeRaster(used_byte, out_used,
    overwrite = TRUE,
    wopt = wopt_byte(opts$SPEED_OVER_SIZE, na = 255L)
  )
}
if (OVERWRITE || !(SKIP_EXISTING && file.exists(out_counts))) {
  writeRaster(c(cnt_cropland, cnt_urban), out_counts,
    overwrite = TRUE,
    wopt = wopt_int(opts$SPEED_OVER_SIZE)
  )
}

tag <- sprintf("glc_used_ge%d_%d-%d", N_YEARS, Y1, Y2)

p_used_global <- tryCatch(global(used_byte, "mean", na.rm = TRUE)[1, 1],
  error = function(e) NA_real_
)

gc()
cat("Done\n")
