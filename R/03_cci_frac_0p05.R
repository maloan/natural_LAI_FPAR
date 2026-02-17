## =============================================================================
# 03_cci_frac_0p05.R
# Aggregate ESA-CCI/C3S land cover to 0.05° fractional cover
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.6)

cci_dir <- cfg$paths$cci_dir
out_dir <- cfg$paths$cci_out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tmpl <- rast(cfg$grids$grid_005$ref_raster)

REMAKE_ALL <- as.logical(Sys.getenv("REMAKE_ALL", "FALSE"))
SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))

ESACCI <- cfg$esa_cci$classes
nodata_vals <- unique(c(ESACCI$nodata, 255))

groups <- Filter(Negate(is.null), list(
  cropland  = ESACCI$cropland,
  urban     = ESACCI$urban,
  cls30     = ESACCI$cls30,
  cls40     = ESACCI$cls40,
  grassland = ESACCI$grassland
))

start_year <- cfg$project$years$cci_start
end_year <- cfg$project$years$cci_end

# Precompute write options once
wopt <- wopt_f32(FALSE)
gdal_opts <- wopt$gdal
naflag <- wopt$NAflag %||% -9999

# ------------------------------------------------------------------------------
# Choose one file per year
# ------------------------------------------------------------------------------
all_files <- list.files(cci_dir, pattern = "\\.tif$", full.names = TRUE)
if (!length(all_files)) {
  stop("No CCI GeoTIFFs found in: ", cci_dir)
}

get_year <- function(x) {
  m <- regexpr("(19|20)\\d{2}", basename(x), perl = TRUE)
  if (m[1] < 0) {
    return(NA_integer_)
  }
  as.integer(substr(basename(x), m[1], m[1] + attr(m, "match.length") - 1))
}
get_source_rank <- function(x) {
  if (grepl("^C3S", basename(x))) {
    2L
  } else {
    1L
  }
}

yrs <- vapply(all_files, get_year, integer(1))
ok <- !is.na(yrs) & yrs >= start_year & yrs <= end_year
all_files <- all_files[ok]
yrs <- yrs[ok]

rank <- vapply(all_files, get_source_rank, integer(1))

# for each year pick file with max rank (C3S preferred); tie-breaker = first
pick_by_year <- tapply(
  seq_along(all_files),
  yrs,
  function(idx) {
    idx[which.max(rank[idx])]
  }
)

plan_year <- as.integer(names(pick_by_year))
plan_path <- all_files[unlist(pick_by_year, use.names = FALSE)]
o_tif <- file.path(out_dir, sprintf("ESACCI_frac_%d_0p05.tif", plan_year))

# ------------------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------------------
for (i in seq_along(plan_year)) {
  yr <- plan_year[i]
  f <- plan_path[i]
  ot <- o_tif[i]

  if (SKIP_EXISTING && file.exists(ot) && !REMAKE_ALL) {
    message("✓ Year ", yr, " already complete — skipping.")
    next
  }

  t0 <- Sys.time()
  message("→ [", yr, "] start")

  r <- rast(f)
  r <- terra::subst(r, nodata_vals, NA)

  m_stack <- rast(lapply(groups, function(cls) {
    classify(r, cbind(cls, 1), others = 0)
  }))

  frac <- resample(m_stack, tmpl, method = "average")
  names(frac) <- paste0("frac_", names(groups))

  w30 <- cfg$esa_cci$weights$cls30 %||% 0.75
  w40 <- cfg$esa_cci$weights$cls40 %||% 0.25

  fc <- frac$frac_cropland %||% 0
  fu <- frac$frac_urban %||% 0
  f30 <- frac$frac_cls30 %||% 0
  f40 <- frac$frac_cls40 %||% 0

  frac_fused <- clamp(fc + fu + w30 * f30 + w40 * f40, 0, 1)
  names(frac_fused) <- "frac_fused"

  frac_grass <- frac$frac_grassland
  names(frac_grass) <- "frac_grass"

  out <- c(frac_fused, frac_grass)

  writeRaster(out, ot, overwrite = TRUE, gdal = gdal_opts, NAflag = naflag)

  rm(r, m_stack, frac, out)
  gc()

  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message("✓ [", yr, "] done (", dt, " s)")
}

message("Done.")
