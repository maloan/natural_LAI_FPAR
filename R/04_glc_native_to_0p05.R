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
test_only <- "--test" %in% commandArgs(trailingOnly = TRUE)

years_5year <- c(1985, 1990, 1995)
years_annual <- 2000:2022

for (command in c("gdal_calc.py", "gdal_translate", "gdalwarp", "gdalbuildvrt", "unzip")) {
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

run_command <- function(command, args, env = character()) {
  status <- system2(command, args, env = env)
  if (status != 0) {
    stop(command, " failed with status ", status)
  }
}

aggregate_tile <- function(file, years, period, test_subset = FALSE) {
  id <- tile_id(file)
  prefix <- if (test_subset) "glc_test" else "glc"
  mode_out <- file.path(tile_out_dir, sprintf("%s_mode_%s_%s_0p05.tif", prefix, period, id))
  grass_out <- file.path(tile_out_dir, sprintf("%s_grass_%s_%s_0p05.tif", prefix, period, id))

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
  nrows_target <- 100
  ncols_target <- 100

  work_dir <- tempfile("glc-tile-")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  if (test_subset) {
    target_extent <- ext(
      xmin(target_extent), xmin(target_extent) + 1,
      ymax(target_extent) - 1, ymax(target_extent)
    )
    source_file <- file.path(work_dir, "source_subset.tif")
    run_command(
      "gdal_translate",
      c(
        "-q", "-of", "GTiff", "-b", "1",
        "-projwin",
        xmin(target_extent), ymax(target_extent),
        xmax(target_extent), ymin(target_extent),
        shQuote(file), shQuote(source_file)
      )
    )
    years <- years[1]
    nrows_target <- 20
    ncols_target <- 20
  }

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
      "-ts", ncols_target, nrows_target,
      "-multi", "-wo", "NUM_THREADS=4",
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

  rm(source_info)
  gc(verbose = FALSE)
}

extract_archive <- function(archive, output_dir, member_pattern = NULL) {
  args <- c("-q", archive)
  if (!is.null(member_pattern)) {
    args <- c(args, member_pattern)
  }
  args <- c(args, "-d", output_dir)

  status <- system2("unzip", args)
  if (status != 0) {
    stop("Could not extract ", basename(archive))
  }
}

process_archive <- function(archive, member_pattern = NULL, test_subset = FALSE) {
  extract_dir <- tempfile("glc-fcs30d-v2-")
  dir.create(extract_dir)
  on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)

  message("→ Extracting ", basename(archive))
  extract_archive(archive, extract_dir, member_pattern)
  files <- list.files(extract_dir, "\\.tif$", recursive = TRUE, full.names = TRUE)

  files_5year <- files[grepl("_5years_", files)]
  files_annual <- files[grepl("_Annual_", files)]
  if (!length(files_5year) || !length(files_annual)) {
    stop("Expected paired 5-year and annual files in ", basename(archive))
  }

  process_files <- function(files, years, period) {
    workers <- if (test_subset) 1L else 3L
    results <- parallel::mclapply(
      files,
      function(file) aggregate_tile(file, years, period, test_subset),
      mc.cores = workers,
      mc.preschedule = FALSE
    )
    failed <- vapply(results, inherits, logical(1), "try-error")
    if (any(failed)) {
      stop("Failed ", period, " tiles:\n", paste(results[failed], collapse = "\n"))
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

  vrt_dir <- tempfile("glc-fcs30d-vrt-")
  dir.create(vrt_dir)
  on.exit(unlink(vrt_dir, recursive = TRUE), add = TRUE)

  build_mosaic <- function(files, output, nodata) {
    input_list <- paste0(output, ".txt")
    writeLines(files, input_list)
    run_command(
      "gdalbuildvrt",
      c(
        "-q", "-overwrite",
        "-resolution", "user", "-tr", "0.05", "0.05",
        "-te", "-180", "-90", "180", "90",
        "-srcnodata", nodata, "-vrtnodata", nodata,
        "-input_file_list", shQuote(input_list),
        shQuote(output)
      )
    )
  }

  mode_5year_vrt <- file.path(vrt_dir, "mode_5year.vrt")
  mode_annual_vrt <- file.path(vrt_dir, "mode_annual.vrt")
  grass_5year_vrt <- file.path(vrt_dir, "grass_5year.vrt")
  grass_annual_vrt <- file.path(vrt_dir, "grass_annual.vrt")

  build_mosaic(mode_5year_files, mode_5year_vrt, "255")
  build_mosaic(mode_annual_files, mode_annual_vrt, "255")
  build_mosaic(grass_5year_files, grass_5year_vrt, "-9999")
  build_mosaic(grass_annual_files, grass_annual_vrt, "-9999")

  years <- c(years_5year, years_annual)
  mode_vrts <- c(rep(mode_5year_vrt, length(years_5year)), rep(mode_annual_vrt, length(years_annual)))
  mode_bands <- c(seq_along(years_5year), seq_along(years_annual))

  for (i in seq_along(years)) {
    year <- years[i]
    output <- file.path(global_out_dir, sprintf("GLC_FCS30D_mode_0p05_%d.tif", year))
    if (!file.exists(output)) {
      run_command(
        "gdal_translate",
        c(
          "-q", "-of", "GTiff", "-b", mode_bands[i],
          "-ot", "Byte", "-a_nodata", "255",
          "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
          shQuote(mode_vrts[i]), shQuote(output)
        )
      )
    }
  }

  grass_output <- file.path(
    global_out_dir,
    "GLC_FCS30D_grass_fraction_0p05_1985_2022.tif"
  )
  if (!file.exists(grass_output)) {
    grass_layer_vrts <- character(length(years))
    grass_vrts <- c(rep(grass_5year_vrt, length(years_5year)), rep(grass_annual_vrt, length(years_annual)))
    grass_bands <- c(seq_along(years_5year), seq_along(years_annual))

    for (i in seq_along(years)) {
      grass_layer_vrts[i] <- file.path(vrt_dir, sprintf("grass_%d.vrt", years[i]))
      run_command(
        "gdal_translate",
        c(
          "-q", "-of", "VRT", "-b", grass_bands[i],
          shQuote(grass_vrts[i]), shQuote(grass_layer_vrts[i])
        )
      )
    }

    grass_list <- file.path(vrt_dir, "grass_layers.txt")
    grass_stack_vrt <- file.path(vrt_dir, "grass_stack.vrt")
    writeLines(grass_layer_vrts, grass_list)
    run_command(
      "gdalbuildvrt",
      c(
        "-q", "-overwrite", "-separate",
        "-srcnodata", "-9999", "-vrtnodata", "-9999",
        "-input_file_list", shQuote(grass_list),
        shQuote(grass_stack_vrt)
      )
    )
    run_command(
      "gdal_translate",
      c(
        "-q", "-of", "GTiff", "-ot", "Float32", "-a_nodata", "-9999",
        "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "TILED=YES",
        shQuote(grass_stack_vrt), shQuote(grass_output)
      )
    )
  }
}

archives <- sort(list.files(
  archive_dir,
  "^GLC_FCS30D_19852022maps_.*\\.zip$",
  full.names = TRUE
))

if (test_only) {
  archive <- file.path(archive_dir, "GLC_FCS30D_19852022maps_E0-E5.zip")
  if (!file.exists(archive)) {
    stop("Representative test archive is missing: ", archive)
  }
  process_archive(archive, "*E0N10*.tif", test_subset = TRUE)
  message("Representative one-year, 1-degree E0N10 aggregation test complete.")
} else {
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
}
