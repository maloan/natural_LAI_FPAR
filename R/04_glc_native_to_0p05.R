# =============================================================================
# 04_glc_native_to_0p05.R — Aggregate native GLC_FCS30D v2 tiles to 0.05°
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "io.R"))

terraOptions(progress = 1, memfrac = 0.25)

archive_dir <- here("data-raw", "GLC_FCS30D", "archives")
tile_out_dir <- here("data", "frac", "glc_native_tiles_0p05")
global_out_dir <- here("data-raw", "GLC_FCS30D")

years_5year <- c(1985, 1990, 1995)
years_annual <- 2000:2022

for (command in c("gdal_calc.py", "gdal_merge.py", "gdal_translate", "gdalwarp", "unzip")) {
  if (!nzchar(Sys.which(command))) {
    stop("Missing required command: ", command)
  }
}

dir.create(tile_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(global_out_dir, recursive = TRUE, showWarnings = FALSE)

tile_extent <- function(file) {
  token <- strcapture(
    "_([EW])([0-9]+)([NS])([0-9]+)_",
    basename(file),
    data.frame(
      lon_dir = character(),
      lon = integer(),
      lat_dir = character(),
      lat = integer()
    )
  )

  if (nrow(token) != 1 || anyNA(token)) {
    stop("Could not read tile coordinates from: ", basename(file))
  }

  west <- if (token$lon_dir == "E") token$lon else -token$lon
  north <- if (token$lat_dir == "N") token$lat else -token$lat

  ext(west, west + 5, north - 5, north)
}

tile_id <- function(file) {
  sub(".*_([EW][0-9]+[NS][0-9]+)_.*", "\\1", basename(file))
}

run_command <- function(command, args) {
  status <- system2(command, args)
  if (status != 0) {
    stop(command, " failed with status ", status)
  }
}

aggregate_tile <- function(file, years, period) {
  id <- tile_id(file)
  mode_out <- file.path(tile_out_dir, sprintf("glc_mode_%s_%s_0p05.tif", period, id))
  grass_out <- file.path(tile_out_dir, sprintf("glc_grass_%s_%s_0p05.tif", period, id))

  if (file.exists(mode_out) && file.exists(grass_out)) {
    return(invisible(NULL))
  }

  message("→ Aggregating ", basename(file))
  source_info <- rast(file)
  if (nlyr(source_info) != length(years)) {
    stop("Unexpected layer count in ", basename(file))
  }

  target_extent <- tile_extent(file)
  source_file <- file
  work_dir <- tempfile("glc-tile-")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  class_file <- file.path(work_dir, "class.tif")
  grass_file <- file.path(work_dir, "grass.tif")

  calculate_raster <- function(expression, output) {
    run_command(
      "gdal_calc.py",
      c(
        "--quiet", "--overwrite", "--hideNoData", "--allBands=A",
        "--type=Byte", "--NoDataValue=255", "--format=GTiff",
        "--co=COMPRESS=DEFLATE", "--co=TILED=YES", "--co=BIGTIFF=IF_SAFER",
        "-A", shQuote(source_file),
        shQuote(paste0("--calc=", expression)),
        "--outfile", shQuote(output)
      )
    )
  }

  calculate_raster("where((A==0)|(A==250),255,A)", class_file)
  calculate_raster("where((A==0)|(A==250),255,A==130)", grass_file)

  warp_args <- function(method, source, output, output_type, nodata) {
    c(
      "-q", "-overwrite", "-r", method,
      "-srcnodata", "255", "-dstnodata", nodata,
      "-te",
      xmin(target_extent), ymin(target_extent),
      xmax(target_extent), ymax(target_extent),
      "-ts", "100", "100",
      "-wo", "NUM_THREADS=1",
      "-ot", output_type,
      "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
      shQuote(source), shQuote(output)
    )
  }
  run_command(
    "gdalwarp",
    warp_args("mode", class_file, mode_out, "Byte", "255")
  )
  run_command(
    "gdalwarp",
    warp_args("average", grass_file, grass_out, "Float32", "-9999")
  )

  # GDAL's average resampling can marginally exceed one at a tile edge.
  # Keep the intermediate fractional layer within its physical range.
  grass_fraction <- clamp(rast(grass_out), lower = 0, upper = 1)
  clamped_file <- file.path(work_dir, "grass_clamped.tif")
  writeRaster(grass_fraction, clamped_file, wopt = wopt_f32(FALSE))
  if (!file.copy(clamped_file, grass_out, overwrite = TRUE)) {
    stop("Could not replace grass-fraction tile: ", grass_out)
  }

  rm(source_info)
  gc(verbose = FALSE)
}

extract_archive <- function(archive, output_dir) {
  status <- system2("unzip", c("-q", archive, "-d", output_dir))
  if (status != 0) {
    stop("Could not extract ", basename(archive))
  }
}

process_archive <- function(archive) {
  extract_dir <- tempfile("glc-fcs30d-v2-")
  dir.create(extract_dir)
  on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)

  message("→ Extracting ", basename(archive))
  extract_archive(archive, extract_dir)
  files <- list.files(extract_dir, "\\.tif$", recursive = TRUE, full.names = TRUE)

  files_5year <- files[grepl("_5years_", files)]
  files_annual <- files[grepl("_Annual_", files)]
  if (!length(files_5year) || !length(files_annual)) {
    stop("Expected paired 5-year and annual files in ", basename(archive))
  }

  process_files <- function(files, years, period) {
    for (file in files) {
      aggregate_tile(file, years, period)
    }
  }

  process_files(files_5year, years_5year, "5year")
  process_files(files_annual, years_annual, "annual")
}

write_global_outputs <- function() {
  mode_5year_files <- sort(list.files(
    tile_out_dir,
    "^glc_mode_5year_.*_0p05\\.tif$",
    full.names = TRUE
  ))
  mode_annual_files <- sort(list.files(
    tile_out_dir,
    "^glc_mode_annual_.*_0p05\\.tif$",
    full.names = TRUE
  ))
  grass_5year_files <- sub("glc_mode_", "glc_grass_", mode_5year_files, fixed = TRUE)
  grass_annual_files <- sub("glc_mode_", "glc_grass_", mode_annual_files, fixed = TRUE)

  if (length(mode_5year_files) != 961 || length(mode_annual_files) != 961) {
    stop(
      "Expected 961 paired native tiles, found ",
      length(mode_5year_files), " five-year and ",
      length(mode_annual_files), " annual tiles."
    )
  }
  stopifnot(all(file.exists(c(grass_5year_files, grass_annual_files))))

  years <- c(years_5year, years_annual)
  mode_outputs <- file.path(
    global_out_dir,
    sprintf("GLC_FCS30D_mode_0p05_%d.tif", years)
  )
  grass_output <- file.path(
    global_out_dir,
    "GLC_FCS30D_grass_fraction_0p05_1985_2022.tif"
  )

  mosaic_dir <- tempfile("glc-fcs30d-mosaic-")
  dir.create(mosaic_dir)
  on.exit(unlink(mosaic_dir, recursive = TRUE), add = TRUE)

  build_mosaic <- function(files, output, nodata, output_type) {
    run_command(
      "gdal_merge.py",
      c(
        "-q", "-o", shQuote(output), "-of", "GTiff",
        "-ul_lr", "-180", "90", "180", "-90",
        "-ps", "0.05", "0.05",
        "-n", nodata, "-init", nodata, "-a_nodata", nodata,
        "-ot", output_type,
        "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
        "-co", "BIGTIFF=YES",
        vapply(files, shQuote, character(1))
      )
    )
  }

  mode_5year_mosaic <- file.path(mosaic_dir, "mode_5year.tif")
  mode_annual_mosaic <- file.path(mosaic_dir, "mode_annual.tif")
  grass_5year_mosaic <- file.path(mosaic_dir, "grass_5year.tif")
  grass_annual_mosaic <- file.path(mosaic_dir, "grass_annual.tif")

  build_mosaic(mode_5year_files, mode_5year_mosaic, "255", "Byte")
  build_mosaic(mode_annual_files, mode_annual_mosaic, "255", "Byte")
  build_mosaic(grass_5year_files, grass_5year_mosaic, "-9999", "Float32")
  build_mosaic(grass_annual_files, grass_annual_mosaic, "-9999", "Float32")

  mode_mosaics <- c(
    rep(mode_5year_mosaic, length(years_5year)),
    rep(mode_annual_mosaic, length(years_annual))
  )
  mode_bands <- c(seq_along(years_5year), seq_along(years_annual))

  for (i in seq_along(years)) {
    output <- mode_outputs[i]
    run_command(
      "gdal_translate",
      c(
        "-q", "-overwrite", "-of", "GTiff", "-b", mode_bands[i],
        "-ot", "Byte", "-a_nodata", "255",
        "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
        shQuote(mode_mosaics[i]), shQuote(output)
      )
    )
  }

  grass <- c(rast(grass_5year_mosaic), rast(grass_annual_mosaic))
  names(grass) <- sprintf("Y%d", years)
  writeRaster(
    grass,
    grass_output,
    overwrite = TRUE,
    wopt = wopt_f32(FALSE)
  )
}

archives <- sort(list.files(
  archive_dir,
  "^GLC_FCS30D_19852022maps_.*\\.zip$",
  full.names = TRUE
))


if (length(archives) != 36) {
    stop("Expected 36 GLC_FCS30D v2 archives in ", archive_dir, "; found ", length(archives))
}

  tile_counts <- vapply(
    c("mode_5year", "grass_5year", "mode_annual", "grass_annual"),
    function(name) {
      length(list.files(tile_out_dir, sprintf("^glc_%s_.*_0p05\\.tif$", name)))
    },
    integer(1)
  )

  if (!all(tile_counts == 961)) {
    for (archive in archives) {
      process_archive(archive)
    }
  }

  write_global_outputs()
  message("GLC_FCS30D v2 local preprocessing complete.")
