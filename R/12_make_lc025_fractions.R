# =============================================================================
# 12_make_lc025_fractions.R — Aggregate annual 300 m ESA-CCI land cover to 0.25°
# =============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(sf)
  library(lwgeom)
  library(Rcpp)
  library(here)
})

source(here("R", "helpers", "landcover_classification.R"))

terraOptions(progress = 1, memfrac = 0.2, todisk = TRUE)

input_dir <- here("data-raw", "ESACCI", "ESACCI_1992-2022")
fraction_dir <- here(
  "analysis", "tmp", paste0("lc025_fraction_yearly_", lc_fraction_schema)
)
years <- 1992:2022

parent_classes <- lc_primary_classes
class_to_parent <- lc_native_to_primary
batch_target_rows <- 2L

input_basename <- function(year) {
  if (year <= 2015) {
    sprintf("ESACCI-LC-L4-LCCS-Map-300m-P1Y-%d-v2.0.7cds.nc", year)
  } else {
    sprintf("C3S-LC-L4-LCCS-Map-300m-P1Y-%d-v2.1.1.nc", year)
  }
}

archive_basename <- function(year) {
  version <- if (year <= 2015) "v2_0_7cds" else "v2_1_1"
  sprintf("C3S_LC_%d_%s.zip", year, version)
}

find_input <- function(year) {
  nc_path <- file.path(input_dir, input_basename(year))
  zip_path <- file.path(input_dir, archive_basename(year))
  if (file.exists(nc_path) && isTRUE(file.info(nc_path)$size > 0)) {
    return(list(type = "netcdf", path = nc_path, member = NA_character_))
  }
  if (!file.exists(zip_path) || !isTRUE(file.info(zip_path)$size > 0)) {
    return(NULL)
  }

  listing <- try(utils::unzip(zip_path, list = TRUE), silent = TRUE)
  if (inherits(listing, "try-error")) {
    return(NULL)
  }
  members <- listing$Name[basename(listing$Name) == input_basename(year)]
  if (length(members) != 1L) {
    return(NULL)
  }
  list(type = "zip", path = zip_path, member = members)
}

validate_netcdf <- function(path) {
  nc <- try(ncdf4::nc_open(path), silent = TRUE)
  if (inherits(nc, "try-error")) {
    return(FALSE)
  }
  on.exit(ncdf4::nc_close(nc))
  valid <- "lccs_class" %in% names(nc$var) &&
    identical(as.integer(nc$var$lccs_class$size), c(129600L, 64800L, 1L))
  isTRUE(valid)
}

# Validate the complete 31-year series before writing derived output.
inputs <- setNames(lapply(years, find_input), years)
missing_inputs <- vapply(inputs, is.null, logical(1))
if (any(missing_inputs)) {
  stop(
    "Missing, incomplete, or invalid ESA-CCI archives/NetCDFs for years: ",
    paste(years[missing_inputs], collapse = ", "),
    call. = FALSE
  )
}
invalid_netcdf <- vapply(inputs, function(input) {
  identical(input$type, "netcdf") && !validate_netcdf(input$path)
}, logical(1))
if (any(invalid_netcdf)) {
  stop(
    "Invalid ESA-CCI NetCDFs for years: ",
    paste(years[invalid_netcdf], collapse = ", "),
    call. = FALSE
  )
}

materialize_input <- function(year) {
  input <- inputs[[as.character(year)]]
  if (identical(input$type, "netcdf")) {
    return(list(path = input$path, temporary_dir = NULL))
  }

  extract_dir <- tempfile(sprintf("esacci_%d_", year))
  dir.create(extract_dir)
  status <- system2(
    "unzip",
    c(
      "-qq", "-o", shQuote(input$path), shQuote(input$member),
      "-d", shQuote(extract_dir)
    )
  )
  extracted <- file.path(extract_dir, input$member)
  if (status != 0L || !file.exists(extracted) || !validate_netcdf(extracted)) {
    unlink(extract_dir, recursive = TRUE)
    stop("Could not extract a valid ESA-CCI NetCDF from ", input$path, call. = FALSE)
  }
  list(path = extracted, temporary_dir = extract_dir)
}

years_to_process <- years[vapply(years, function(year) {
  fraction_file <- file.path(fraction_dir, sprintf("lc025_fraction_%d.tif", year))
  !lc_fraction_file_is_valid(fraction_file) ||
    file.mtime(fraction_file) < file.mtime(inputs[[as.character(year)]]$path)
}, logical(1))]
if (!length(years_to_process)) {
  message("All annual land-cover fractions are current.")
  quit(save = "no", status = 0L)
}

lookup <- integer(256)
class_index <- setNames(seq_along(parent_classes), parent_classes)
lookup[as.integer(names(class_to_parent)) + 1L] <- unname(class_index[as.character(class_to_parent)])

# The explicit polygons reproduce the WGS84 row-area weights used previously.
latitude_edges <- 90 - (0:64800) / 360
row_polygons <- lapply(seq_len(64800), function(row) {
  sf::st_polygon(list(matrix(
    c(
      0, latitude_edges[row],
      1 / 360, latitude_edges[row],
      1 / 360, latitude_edges[row + 1L],
      0, latitude_edges[row + 1L],
      0, latitude_edges[row]
    ),
    ncol = 2,
    byrow = TRUE
  )))
})
row_area <- as.numeric(lwgeom::st_geod_area(sf::st_sfc(row_polygons, crs = 4326)))
rm(row_polygons)

Rcpp::cppFunction('
NumericMatrix aggregate_land_cover(
    IntegerMatrix source,
    IntegerVector lookup,
    NumericVector row_area,
    int class_count) {
  int source_columns = source.nrow();
  int source_rows = source.ncol();
  int target_columns = source_columns / 90;
  NumericMatrix fractions(class_count, target_columns);
  NumericVector totals(target_columns);

  for (int row = 0; row < source_rows; ++row) {
    double weight = row_area[row];
    for (int column = 0; column < source_columns; ++column) {
      int value = source(column, row);
      if (value < 0 || value >= lookup.size()) continue;
      int class_id = lookup[value];
      if (class_id == 0) continue;
      int target_column = column / 90;
      fractions(class_id - 1, target_column) += weight;
      totals[target_column] += weight;
    }
  }

  for (int column = 0; column < target_columns; ++column) {
    if (totals[column] == 0) {
      for (int class_id = 0; class_id < class_count; ++class_id) {
        fractions(class_id, column) = NA_REAL;
      }
    } else {
      for (int class_id = 0; class_id < class_count; ++class_id) {
        fractions(class_id, column) /= totals[column];
      }
    }
  }
  return fractions;
}
')

aggregate_year <- function(year) {
  materialized <- materialize_input(year)
  path <- materialized$path
  on.exit({
    if (!is.null(materialized$temporary_dir)) {
      unlink(materialized$temporary_dir, recursive = TRUE)
    }
  }, add = TRUE)

  fractions <- array(
    NA_real_,
    dim = c(length(parent_classes), 720L, 1440L)
  )

  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expected_size <- c(129600L, 64800L, 1L)
  if (!identical(nc$var$lccs_class$size, expected_size)) {
    stop("Unexpected lccs_class dimensions in ", path)
  }

  # Read adjacent rows together to avoid repeated NetCDF decompression.
  for (batch_start in seq.int(1L, 720L, by = batch_target_rows)) {
    batch_end <- min(batch_start + batch_target_rows - 1L, 720L)
    source_start <- (batch_start - 1L) * 90L + 1L
    source_batch <- ncdf4::ncvar_get(
      nc,
      "lccs_class",
      start = c(1L, source_start, 1L),
      count = c(129600L, (batch_end - batch_start + 1L) * 90L, 1L),
      collapse_degen = TRUE,
      raw_datavals = TRUE
    )
    for (target_row in batch_start:batch_end) {
      offset <- (target_row - batch_start) * 90L
      source <- source_batch[, (offset + 1L):(offset + 90L), drop = FALSE]
      row_start <- (target_row - 1L) * 90L + 1L
      fractions[, target_row, ] <- aggregate_land_cover(
        source,
        lookup,
        row_area[row_start:(row_start + 89L)],
        length(parent_classes)
      )

      if (target_row %% 20L == 0L || target_row == 720L) {
        message("  target rows ", target_row - 19L, "-", target_row)
      }
    }
    rm(source_batch, source)
  }
  fractions
}

write_outputs <- function(year, fractions) {
  fraction_values <- vapply(
    seq_along(parent_classes),
    function(class_id) as.vector(t(fractions[class_id, , ])),
    numeric(720L * 1440L)
  )

  fraction_raster <- rast(
    nrows = 720,
    ncols = 1440,
    nlyrs = length(parent_classes),
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 90,
    crs = "EPSG:4326"
  )
  values(fraction_raster) <- fraction_values
  names(fraction_raster) <- paste0("lc_", parent_classes)

  output_file <- file.path(fraction_dir, sprintf("lc025_fraction_%d.tif", year))
  tmp <- paste0(output_file, ".tmp.tif")
  writeRaster(
    fraction_raster,
    tmp,
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
  )
  if (!file.rename(tmp, output_file)) {
    stop("Could not replace ", output_file)
  }
}

dir.create(fraction_dir, recursive = TRUE, showWarnings = FALSE)

for (year in years_to_process) {
  message(
    "Processing ", year, ": ",
    basename(inputs[[as.character(year)]]$path)
  )
  fractions <- aggregate_year(year)
  write_outputs(year, fractions)
  rm(fractions)
  gc(verbose = FALSE)
  message("Completed ", year)
}
