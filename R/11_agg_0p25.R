## =============================================================================
# 11_agg_0p25.R — Area-weighted aggregation of masked LAI/FPAR (0.05° → 0.25°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

# --- repo helpers --------------------------------------------------------------
source(here("R", "helpers", "paths.R"))
source(here("R", "helpers", "files.R"))
source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))
source(here("R", "helpers", "plotting.R"))

cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

# ------------------------------------------------------------------------------
# Env
# ------------------------------------------------------------------------------
skip_existing <- as_bool(Sys.getenv("skip_existing"), default = FALSE)
overwrite <- as_bool(Sys.getenv("overwrite"), default = FALSE)
remake_ql <- as_bool(Sys.getenv("remake_ql"), default = FALSE)

var <- toupper(Sys.getenv("var", "FPAR")) # LAI|FPAR
mask <- toupper(Sys.getenv("mask", "CCI")) # CCI|GLC
if (!var %in% c("LAI", "FPAR")) {
  stop("Unsupported var: ", var, ". Use LAI or FPAR")
}
if (!mask %in% c("CCI", "GLC")) {
  stop("Unsupported mask: ", mask, ". Use CCI or GLC")
}

# ------------------------------------------------------------------------------
# Refs and weights
# ------------------------------------------------------------------------------
ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)

if (!compareGeom(area005, ref005, stopOnError = FALSE)) {
  area005 <- resample(area005, ref005, method = "bilinear")
}

# ------------------------------------------------------------------------------
# dirs
# ------------------------------------------------------------------------------
key_in <- sprintf("masked_%s_%s_0p05_dir", tolower(var), tolower(mask))
key_out <- sprintf("masked_%s_%s_0p25_dir", tolower(var), tolower(mask))

in_dir <- cfg$paths[[key_in]]
out_dir <- cfg$paths[[key_out]]

stopifnot(is.character(in_dir), length(in_dir) == 1, nzchar(in_dir))
stopifnot(is.character(out_dir), length(out_dir) == 1, nzchar(out_dir))

patt <- sprintf("^%s_\\d{6}_0p05_masked\\.tif$", var)
ql_title <- var

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
qdir <- file.path(out_dir, "quicklooks")
dir.create(qdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# inputs
# ------------------------------------------------------------------------------
stopifnot(dir.exists(in_dir))
files <- sort(list.files(in_dir, pattern = patt, full.names = TRUE))
stopifnot(length(files) > 0L)

# ------------------------------------------------------------------------------
# Aggregation loop
# ------------------------------------------------------------------------------
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_masked_%s_0p25.tif", var, ym))
  do_write <- overwrite || !file.exists(out) || !skip_existing
  do_ql <- (substr(ym, 5, 6) %in% c("01", "07")) &&
    (remake_ql || !file.exists(file.path(
      qdir, sprintf("quicklook_%s_0p25_%s.png", ql_title, ym)
    )))

  if (!do_write && !do_ql) {
    next
  }

  if (do_write) {
    r <- rast(f)

    if (!compareGeom(r, ref005, stopOnError = FALSE)) {
      r <- resample(r, ref005, method = "bilinear")
    }

    # numerator/denominator (area-weighted mean, handling NA)
    num <- aggregate(r * area005, fact = 5, fun = "sum", na.rm = TRUE)
    den <- aggregate((!is.na(r)) * area005, fact = 5, fun = "sum", na.rm = TRUE)

    r025 <- ifel(den == 0, NA, num / den)

    if (!compareGeom(r025, ref025, stopOnError = FALSE)) {
      r025 <- resample(r025, ref025, method = "near")
    }
    wopt <- wopt_f32(opts$speed_over_size)
    writeRaster(
      r025,
      out,
      overwrite = TRUE,
      wopt = wopt
    )
  } else {
    r025 <- rast(out)
    if (!compareGeom(r025, ref025, stopOnError = FALSE)) {
      r025 <- resample(r025, ref025, method = "near")
    }
  }

  if (do_ql) {
    out_png <- file.path(qdir, sprintf("quicklook_%s_0p25_%s.png", ql_title, ym))
    write_quicklook_raster(
      r = r025,
      out_png = out_png,
      title = sprintf("%s 0.25° %s", ql_title, ym),
      width = 1400,
      height = 700,
      res = 120,
      legend = TRUE
    )
  }

  gc()
}

gc()
message("Aggregated monthly ", var, " written to: ", out_dir)
