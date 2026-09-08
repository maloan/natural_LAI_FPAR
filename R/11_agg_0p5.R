# =============================================================================
# 11_agg_0p5.R — Area-weighted aggregation of masked LAI/FPAR (0.05° → 0.5°)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})


source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.25)


var <- toupper(Sys.getenv("VAR", "LAI")) # allowed values: LAI or FPAR
mask <- toupper(Sys.getenv("MASK", "CCI")) # allowed values: CCI or GLC


# Refs and weights
ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref050 <- aggregate(ref005,
  fact = 10,
  fun = "mean",
  na.rm = TRUE
) # template
area005 <- rast(cfg$grids$grid_005$area_raster)


# dirs
key_in <- sprintf("masked_%s_%s_0p05_dir", tolower(var), tolower(mask))
key_out <- sprintf("masked_%s_%s_0p5_dir", tolower(var), tolower(mask))

in_dir <- cfg$paths[[key_in]]
out_dir <- cfg$paths[[key_out]]

stopifnot(
  is.character(in_dir),
  length(in_dir) == 1,
  nzchar(in_dir),
  dir.exists(in_dir)
)
stopifnot(is.character(out_dir), length(out_dir) == 1, nzchar(out_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# inputs
stopifnot(dir.exists(in_dir))
months <- format(
  seq(
    as.Date(sprintf("%d-01-01", cfg$project$years$lai_start)),
    as.Date(sprintf("%d-12-01", cfg$project$years$lai_end)),
    by = "month"
  ),
  "%Y%m"
)
files <- file.path(in_dir, sprintf("%s_%s_0p05_masked.tif", var, months))
missing_files <- files[!file.exists(files)]
if (length(missing_files)) {
  stop("Missing expected monthly masked files:\n", paste(missing_files, collapse = "\n"))
}

#  loop
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_masked_%s_0p5.tif", var, ym))
  do_write <- !file.exists(out)
  if (!do_write) {
    next
  }

  r <- rast(f)
  r <- align_to_template(r, ref005, method = "bilinear")

  num <- aggregate(r * area005,
    fact = 10,
    fun = "sum",
    na.rm = TRUE
  )
  den <- aggregate((!is.na(r)) * area005,
    fact = 10,
    fun = "sum",
    na.rm = TRUE
  )
  r050 <- ifel(den == 0, NA, num / den)
  r050 <- align_to_template(r050, ref050, method = "near")

  wopt <- wopt_f32(FALSE)
  writeRaster(r050, out, overwrite = TRUE, wopt = wopt)
}
