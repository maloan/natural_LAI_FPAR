# =============================================================================
# 11_agg_0p25.R — Area-weighted aggregation of masked LAI/FPAR (0.05° → 0.25°)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})


source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "plotting.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.25)


var <- toupper(Sys.getenv("VAR", "LAI")) # allowed values: LAI or FPAR
mask <- toupper(Sys.getenv("MASK", "CCI")) # allowed values: CCI or GLC
if (!var %in% c("LAI", "FPAR")) {
  stop_msg("Unsupported var: ", var, ". Use LAI or FPAR")
}
if (!mask %in% c("CCI", "GLC")) {
  stop_msg("Unsupported mask: ", mask, ". Use CCI or GLC")
}


# Refs and weights
ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)

area005 <- align_to_template(area005, ref005, method = "bilinear")

# dirs
key_in <- sprintf("masked_%s_%s_0p05_dir", tolower(var), tolower(mask))
key_out <- sprintf("masked_%s_%s_0p25_dir", tolower(var), tolower(mask))

in_dir <- cfg$paths[[key_in]]
out_dir <- cfg$paths[[key_out]]

stopifnot(is.character(in_dir), length(in_dir) == 1, nzchar(in_dir))
stopifnot(is.character(out_dir), length(out_dir) == 1, nzchar(out_dir))

ql_title <- var

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
qdir <- file.path(out_dir, "quicklooks")
dir.create(qdir, recursive = TRUE, showWarnings = FALSE)


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

# Aggregation loop
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_masked_%s_0p25.tif", var, ym))

  do_write <- !file.exists(out)
  do_ql <- (substr(ym, 5, 6) %in% c("01", "07")) &&
    (!file.exists(file.path(
      qdir, sprintf("quicklook_%s_0p25_%s.png", ql_title, ym)
    )))

  if (!do_write && !do_ql) {
    next
  }

  if (do_write) {
    r <- rast(f)
    r <- align_to_template(r, ref005, method = "bilinear")

    # numerator/denominator (area-weighted mean, handling NA)
    num <- aggregate(r * area005,
      fact = 5,
      fun = "sum",
      na.rm = TRUE
    )
    den <- aggregate((!is.na(r)) * area005,
      fact = 5,
      fun = "sum",
      na.rm = TRUE
    )
    r025 <- ifel(den == 0, NA, num / den)
    r025 <- align_to_template(r025, ref025, method = "near")

    wopt <- wopt_f32(FALSE)
    writeRaster(r025, out, overwrite = TRUE, wopt = wopt)
  } else {
    # Load existing aggregated file for quicklook
    r025 <- rast(out)
    r025 <- align_to_template(r025, ref025, method = "near")
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

  rm(r025)
  if (do_write) {
    rm(r, num, den)
  }
  gc(verbose = FALSE)
}

gc()
