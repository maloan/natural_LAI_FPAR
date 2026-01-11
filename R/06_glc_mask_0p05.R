## =============================================================================
# 06_glc_mask_0p05.R — Build “used ≥ N years” mask from GLC_FCS30D yearstack
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
# source(here("R","helpers", "geom.R"))
source(here("R", "helpers", "viz.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.9)

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
OVERWRITE <- as.logical(Sys.getenv("OVERWRITE", "FALSE"))
REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))

N_YEARS <- as.integer(Sys.getenv("USED_N_YEARS", "3"))
if (!is.finite(N_YEARS) || N_YEARS < 1) {
  stop("USED_N_YEARS must be a positive integer")
}

glc_out_dir <- cfg$paths$glc_out_dir
masks_dir <- cfg$paths$masks_glc_dir
ql_dir <- file.path(masks_dir, "quicklooks")
dir.create(masks_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

tmpl <- rast(cfg$grids$grid_005$ref_raster)

classes <- cfg$glc$classes
vec_int <- function(x) as.integer(unlist(x, use.names = FALSE))
cropland_vals <- vec_int(classes$cropland)
urban_vals <- vec_int(classes$urban)
nodata_vals <- vec_int(classes$nodata)

# analysis window (match CCI window)
yrs <- cfg$project$years$cci_start:cfg$project$years$cci_end
Y1 <- min(yrs)
Y2 <- max(yrs)

stack_path <- file.path(glc_out_dir, "glc_cat_yearstack_0p05.tif")
if (!file.exists(stack_path)) stop("GLC yearstack not found: ", stack_path)
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
ql_probe <- file.path(
  ql_dir, "global",
  sprintf("quicklook_mask_global_%s.png", tag)
)
if (REMAKE_QL || !file.exists(ql_probe)) {
  quicklook_mask_all_aois(
    mask = used_byte,
    title = sprintf("GLC mask (used ≥ %d years)", N_YEARS),
    tag = tag,
    cfg = cfg,
    ql_root = ql_dir,
    down = 4L,
    include_global = TRUE,
    drop_global_key = FALSE
  )
}

p_used_global <- tryCatch(global(used_byte, "mean", na.rm = TRUE)[1, 1],
  error = function(e) NA_real_
)

cat(sprintf(
  "Wrote:\n  - %s\n  - %s\nQuicklooks in: %s\nUsed≥%d proportion (global):
  %.4f\n",
  out_used, out_counts, ql_dir, N_YEARS, p_used_global
))

gc()
cat("Done\n")
