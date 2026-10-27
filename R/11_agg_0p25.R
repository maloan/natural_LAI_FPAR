# =============================================================================
# 11_agg_0p25.R — Area-weighted aggregation of masked LAI/FPAR (0.05° → 0.25°)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})


source(here("R", "helpers", "io.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.25)


var <- toupper(Sys.getenv("VAR")) # allowed values: LAI or FPAR
mask <- if (cfg$project$run_tag == "glc") "GLC" else "CCI"
if (!var %in% c("LAI", "FPAR")) {
  stop_msg("Unsupported var: ", var, ". Use LAI or FPAR")
}
scenario_output_root(mask, cfg$project$run_tag)


# Refs and weights
ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)

compareGeom(area005, ref005, stopOnError = TRUE)

# dirs
key_in <- sprintf("masked_%s_%s_0p05_dir", tolower(var), tolower(mask))
key_out <- sprintf("masked_%s_%s_0p25_dir", tolower(var), tolower(mask))

in_dir <- cfg$paths[[key_in]]
out_dir <- cfg$paths[[key_out]]

stopifnot(is.character(in_dir), length(in_dir) == 1, nzchar(in_dir))
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
source_dependencies <- c(
  cfg$grids$grid_005$area_raster,
  cfg$grids$grid_005$ref_raster,
  cfg$grids$grid_025$ref_raster,
  here("config", sprintf("config_%s.yml", cfg$project$run_tag)),
  here("R", "11_agg_0p25.R"),
  here("R", "helpers", "io.R")
)
missing_dependencies <- source_dependencies[!file.exists(source_dependencies)]
if (length(missing_dependencies)) {
  stop("Missing aggregation dependencies:\n",
    paste(missing_dependencies, collapse = "\n"))
}
source_mtime <- max(file.mtime(source_dependencies))
wopt <- wopt_f32(FALSE)

# Aggregation loop
for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_masked_%s_0p25.tif", var, ym))

  if (file.exists(out) && isTRUE(file.info(out)$size > 0) &&
      file.mtime(out) >= max(file.mtime(f), source_mtime) &&
      isTRUE(tryCatch(compareGeom(rast(out), ref025,
        stopOnError = FALSE), error = function(e) FALSE))) {
    next
  }

  r <- rast(f)
  compareGeom(r, ref005, stopOnError = TRUE)

  # Mean over valid fine cells, weighted by their surface area.
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
  compareGeom(r025, ref025, stopOnError = TRUE)

  tmp <- paste0(out, ".tmp.tif")
  writeRaster(r025, tmp, overwrite = TRUE, wopt = wopt)
  if (!file.rename(tmp, out)) {
    stop("Could not replace ", out)
  }

  rm(r, num, den, r025)
  gc(verbose = FALSE)
}

gc()
