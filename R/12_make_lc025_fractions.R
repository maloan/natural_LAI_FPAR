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

terraOptions(progress = 1, memfrac = 0.6)

input_dir <- here("data-raw", "ESACCI", "ESACCI_1992-2022")
fraction_dir <- here("analysis", "tmp", "lc025_fraction_yearly")
years <- 1992:2022

parent_classes <- c(
  10L, 50L, 60L, 70L, 80L, 90L, 100L, 120L, 130L,
  140L, 150L, 160L, 180L, 190L, 200L, 210L, 220L
)

class_to_parent <- c(
  `10` = 10L, `11` = 10L, `12` = 10L, `20` = 10L,
  `30` = 10L, `40` = 10L, `50` = 50L,
  `60` = 60L, `61` = 60L, `62` = 60L,
  `70` = 70L, `71` = 70L, `72` = 70L,
  `80` = 80L, `81` = 80L, `82` = 80L,
  `90` = 90L, `100` = 100L, `110` = 100L,
  `120` = 120L, `121` = 120L, `122` = 120L,
  `130` = 130L, `140` = 140L,
  `150` = 150L, `151` = 150L, `152` = 150L, `153` = 150L,
  `160` = 160L, `170` = 160L, `180` = 180L,
  `190` = 190L, `200` = 200L, `201` = 200L, `202` = 200L,
  `210` = 210L, `220` = 220L
)

input_file <- function(year) {
  filename <- if (year <= 2015) {
    sprintf("ESACCI-LC-L4-LCCS-Map-300m-P1Y-%d-v2.0.7cds.nc", year)
  } else {
    sprintf("C3S-LC-L4-LCCS-Map-300m-P1Y-%d-v2.1.1.nc", year)
  }
  file.path(input_dir, filename)
}

missing_inputs <- vapply(years, function(year) !file.exists(input_file(year)), logical(1))
if (any(missing_inputs)) {
  stop("Missing expected ESA-CCI files for years: ", paste(years[missing_inputs], collapse = ", "))
}
years_to_process <- years[vapply(years, function(year) {
  fraction_file <- file.path(fraction_dir, sprintf("lc025_fraction_%d.tif", year))
  !file.exists(fraction_file) || file.mtime(fraction_file) < file.mtime(input_file(year))
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
  path <- input_file(year)
  if (!file.exists(path)) {
    stop("Missing expected ESA-CCI file: ", path)
  }

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

  for (target_row in seq_len(720L)) {
    source_start <- (target_row - 1L) * 90L + 1L
    source <- ncdf4::ncvar_get(
      nc,
      "lccs_class",
      start = c(1L, source_start, 1L),
      count = c(129600L, 90L, 1L),
      collapse_degen = TRUE,
      raw_datavals = TRUE
    )
    source <- matrix(as.integer(source), nrow = 129600L, ncol = 90L)
    fractions[, target_row, ] <- aggregate_land_cover(
      source,
      lookup,
      row_area[source_start:(source_start + 89L)],
      length(parent_classes)
    )

    if (target_row %% 20L == 0L || target_row == 720L) {
      message("  target rows ", target_row - 19L, "-", target_row)
    }
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

  writeRaster(
    fraction_raster,
    file.path(fraction_dir, sprintf("lc025_fraction_%d.tif", year)),
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
  )
}

dir.create(fraction_dir, recursive = TRUE, showWarnings = FALSE)

for (year in years_to_process) {
  message("Processing ", year, ": ", basename(input_file(year)))
  fractions <- aggregate_year(year)
  write_outputs(year, fractions)
  rm(fractions)
  gc(verbose = FALSE)
  message("Completed ", year)
}
