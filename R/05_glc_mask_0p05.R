# =============================================================================
# 05_glc_mask_0p05.R — Build “used ≥ N years” mask from GLC_FCS30D yearstack
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "paths.R"))
source(here("R", "helpers", "files.R"))
source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "plotting.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.9)

skip_existing <- as_bool(Sys.getenv("skip_existing"), default = TRUE)
overwrite <- as_bool(Sys.getenv("overwrite"), default = FALSE)
n_years <- as.integer(Sys.getenv("used_n_years", "3"))

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
year_1 <- min(yrs)
year_2 <- max(yrs)

stack_path <- file.path(glc_out_dir, "glc_cat_yearstack_0p05.tif")
if (!file.exists(stack_path)) {
  stop("GLC yearstack not found: ", stack_path)
}
message("Using GLC yearstack: ", stack_path)

s <- rast(stack_path)
if (is.na(crs(s))) {
  crs(s) <- crs(tmpl)
}
s <- align_to_template(s, tmpl, method = "near")
if (length(nodata_vals)) {
  s[s %in% nodata_vals] <- NA
}

# Extract years from layer names (format: "Y1990")
layer_years <- as.integer(substr(names(s), 2, 5))
keep_idx <- which(layer_years %in% yrs)
s <- s[[keep_idx]]
stopifnot(nlyr(s) > 0)

message(sprintf(
  "Year window for used-counts: %d..%d (%d years)",
  year_1,
  year_2,
  nlyr(s)
))

ncores <- opts$n_workers

cnt_cropland <- app(s, function(v, vals) {
  sum(v %in% vals, na.rm = TRUE)
}, vals = cropland_vals, cores = ncores)
cnt_urban <- app(s, function(v, vals) {
  sum(v %in% vals, na.rm = TRUE)
}, vals = urban_vals, cores = ncores)

names(cnt_cropland) <- "cnt_cropland"
names(cnt_urban) <- "cnt_urban"

cnt_total <- cnt_cropland + cnt_urban
used_byte <- ifel(cnt_total >= n_years, 1L, 0L)

out_used <- file.path(
  masks_dir,
  sprintf("mask_used_ge%d_%d-%d_0p05.tif", n_years, year_1, year_2)
)
out_counts <- file.path(masks_dir, "glc_counts_crop_urban_0p05.tif")

if (overwrite || !(skip_existing && file.exists(out_used))) {
  writeRaster(
    used_byte,
    out_used,
    overwrite = TRUE,
    wopt = wopt_byte(opts$speed_over_size, na = 255L)
  )
}
if (overwrite || !(skip_existing && file.exists(out_counts))) {
  writeRaster(
    c(cnt_cropland, cnt_urban),
    out_counts,
    overwrite = TRUE,
    wopt = wopt_int(opts$speed_over_size)
  )
}

gc()
