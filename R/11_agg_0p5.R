# =============================================================================
# 11_agg_0p5.R — Masked CCI alpha 0.1 fAPAR at 0.5°
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
  library(ncdf4)
})

source(here("R", "helpers", "io.R"))

cfg <- cfg_read()
if (cfg$project$run_tag != "alpha_0.1") {
  stop("This product uses the CCI alpha 0.1 and pasture mask; set RUN_TAG=alpha_0.1.")
}

terraOptions(progress = 1, memfrac = 0.2)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)
compareGeom(area005, ref005, stopOnError = TRUE)
if (nrow(ref005) %% 10L != 0L || ncol(ref005) %% 10L != 0L ||
    !isTRUE(all.equal(res(ref005), c(0.05, 0.05)))) {
  stop("The 0.05-degree reference grid cannot be aggregated by a factor of 10.")
}

ref05 <- rast(
  nrows = nrow(ref005) / 10L,
  ncols = ncol(ref005) / 10L,
  xmin = xmin(ref005), xmax = xmax(ref005),
  ymin = ymin(ref005), ymax = ymax(ref005),
  crs = crs(ref005)
)

in_dir <- cfg$paths$masked_fpar_cci_0p05_dir
if (!is.character(in_dir) || length(in_dir) != 1L || !dir.exists(in_dir)) {
  stop("Missing masked 0.05-degree CCI fAPAR input directory.")
}
out_dir <- here(
  "output", cfg$project$run_tag, "masked_0p5", "FPAR", "masked_FPAR_CCI"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
netcdf_dir <- here("output", "chapter2", "fpar_masked")
dir.create(netcdf_dir, recursive = TRUE, showWarnings = FALSE)

months <- format(
  seq(
    as.Date(sprintf("%d-01-01", cfg$project$years$lai_start)),
    as.Date(sprintf("%d-12-01", cfg$project$years$lai_end)),
    by = "month"
  ),
  "%Y%m"
)
input_files <- file.path(in_dir, sprintf("FPAR_%s_0p05_masked.tif", months))
monthly_files <- file.path(out_dir, sprintf("FPAR_masked_%s_0p5.tif", months))
netcdf_file <- file.path(
  netcdf_dir,
  sprintf("fpar_masked_0p5_monthly_%d-%d.nc",
    cfg$project$years$lai_start, cfg$project$years$lai_end)
)

missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop("Missing masked fAPAR inputs:\n", paste(missing_inputs, collapse = "\n"))
}
source_dependencies <- c(
  cfg$grids$grid_005$ref_raster,
  cfg$grids$grid_005$area_raster,
  here("config", "config_alpha_0.1.yml"),
  here("R", "11_agg_0p5.R"),
  here("R", "helpers", "io.R")
)
missing_dependencies <- source_dependencies[!file.exists(source_dependencies)]
if (length(missing_dependencies)) {
  stop("Missing aggregation dependencies:\n",
    paste(missing_dependencies, collapse = "\n"))
}
source_mtime <- max(file.mtime(source_dependencies))
wopt <- wopt_f32(FALSE)

for (i in seq_along(months)) {
  out <- monthly_files[i]
  if (file.exists(out) && isTRUE(file.info(out)$size > 0) &&
      file.mtime(out) >= max(file.mtime(input_files[i]), source_mtime) &&
      isTRUE(tryCatch(compareGeom(rast(out), ref05,
        stopOnError = FALSE), error = function(e) FALSE))) {
    next
  }

  message("Aggregating masked fAPAR ", months[i], " (", i, "/", length(months), ")")
  r <- rast(input_files[i])
  compareGeom(r, ref005, stopOnError = TRUE)
  if (nlyr(r) != 1L) {
    stop("Expected one fAPAR layer in ", input_files[i])
  }

  # Mean over valid fine cells, weighted by their surface area.
  num <- aggregate(r * area005, fact = 10, fun = "sum", na.rm = TRUE)
  den <- aggregate((!is.na(r)) * area005,
    fact = 10, fun = "sum", na.rm = TRUE)
  r05 <- ifel(den > 0, num / den, NA)
  compareGeom(r05, ref05, stopOnError = TRUE)

  tmp <- paste0(out, ".tmp.tif")
  writeRaster(r05, tmp, overwrite = TRUE, wopt = wopt)
  if (!file.rename(tmp, out)) {
    stop("Could not replace ", out)
  }
  rm(r, num, den, r05)
  gc(verbose = FALSE)
}

if (!file.exists(netcdf_file) ||
    file.mtime(netcdf_file) < max(file.mtime(monthly_files), source_mtime)) {
  message("Writing 516-month masked fAPAR NetCDF.")
  dates <- as.Date(paste0(months, "01"), format = "%Y%m%d")
  lon <- ncdim_def("lon", "degrees_east",
    vals = xFromCol(ref05, seq_len(ncol(ref05))))
  lat <- ncdim_def("lat", "degrees_north",
    vals = yFromRow(ref05, seq_len(nrow(ref05))))
  month_time <- ncdim_def("time", "days since 1970-01-01",
    vals = as.numeric(dates - as.Date("1970-01-01")))
  fpar_var <- ncvar_def(
    "fpar", "1", list(lon, lat, month_time),
    missval = -9999,
    longname = "fAPAR after CCI alpha 0.1, pasture, and nonvegetated exclusions",
    prec = "float", compression = 4,
    chunksizes = c(ncol(ref05), nrow(ref05), 1L)
  )
  tmp_nc <- paste0(netcdf_file, ".tmp.nc")
  nc <- nc_create(tmp_nc, list(fpar_var), force_v4 = TRUE)
  ncatt_put(nc, 0, "Conventions", "CF-1.8")
  ncatt_put(nc, "lon", "standard_name", "longitude")
  ncatt_put(nc, "lat", "standard_name", "latitude")
  ncatt_put(nc, "time", "standard_name", "time")
  for (i in seq_along(monthly_files)) {
    r05 <- rast(monthly_files[i])
    compareGeom(r05, ref05, stopOnError = TRUE)
    vals <- values(r05)[, 1]
    vals[is.na(vals)] <- -9999
    ncvar_put(nc, "fpar",
      array(vals, dim = c(ncol(ref05), nrow(ref05), 1L)),
      start = c(1L, 1L, i),
      count = c(ncol(ref05), nrow(ref05), 1L))
  }
  nc_close(nc)
  metadata <- nc_open(tmp_nc)
  if (!all(c("lon", "lat", "time") %in% names(metadata$dim)) ||
      !"fpar" %in% names(metadata$var)) {
    nc_close(metadata)
    stop("The masked fAPAR NetCDF has incorrect variable or dimension names.")
  }
  nc_close(metadata)
  check <- rast(tmp_nc)
  if (nlyr(check) != length(months) ||
      !identical(as.Date(time(check)), dates)) {
    stop("The masked fAPAR NetCDF has incorrect monthly dates.")
  }
  compareGeom(check[[1]], ref05, stopOnError = TRUE)
  rm(check)
  if (!file.rename(tmp_nc, netcdf_file)) {
    stop("Could not replace ", netcdf_file)
  }
}

message("Masked 0.5-degree fAPAR monthly rasters: ", out_dir)
message("Chapter 2 masked fAPAR NetCDF: ", netcdf_file)
