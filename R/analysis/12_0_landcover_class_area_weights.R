# ==============================================================================
# 12_0_landcover_class_area_weights.R — Exact 0.25-degree land-cover class
# areas before and after masking, derived on the 0.05-degree mask grid
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(Rcpp)
  library(terra)
  library(readr)
  library(tibble)
  library(here)
})

source(here("R", "helpers", "cli_args.R"))
source(here("R", "helpers", "io.R"))

terraOptions(progress = 1, memfrac = 0.45, todisk = TRUE)
if (!nzchar(Sys.getenv("OMP_NUM_THREADS"))) {
  detected_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(detected_cores)) detected_cores <- 1L
  Sys.setenv(OMP_NUM_THREADS = min(8L, detected_cores))
}

default_cfg <- list(
  mask = "CCI",
  alpha = "alpha_0.1",
  lc_year_start = 1992L,
  lc_year_end = 2022L
)
cfg <- parse_cli_args(default_cfg)

mask_source <- toupper(as.character(cfg$mask))
run_tag <- as.character(cfg$alpha)
years <- as.integer(cfg$lc_year_start):as.integer(cfg$lc_year_end)

if (!mask_source %in% c("CCI", "GLC")) {
  stop("mask must be CCI or GLC", call. = FALSE)
}
if (length(years) < 1L || any(!is.finite(years))) {
  stop("Invalid land-cover year range", call. = FALSE)
}

input_dir <- here("data-raw", "ESACCI", "ESACCI_1992-2022")
out_dir <- here("analysis", "tmp", "lc025_class_area_weights")
qa_dir <- here("analysis", "results", "tables", "land_cover")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

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

lookup <- integer(256L)
class_index <- setNames(seq_along(parent_classes), parent_classes)
lookup[as.integer(names(class_to_parent)) + 1L] <-
  unname(class_index[as.character(class_to_parent)])

input_file <- function(year) {
  filename <- if (year <= 2015L) {
    sprintf("ESACCI-LC-L4-LCCS-Map-300m-P1Y-%d-v2.0.7cds.nc", year)
  } else {
    sprintf("C3S-LC-L4-LCCS-Map-300m-P1Y-%d-v2.1.1.nc", year)
  }
  file.path(input_dir, filename)
}

scenario_token <- paste0(tolower(mask_source), "_", gsub("\\.", "p", run_tag))
unmasked_path <- file.path(
  out_dir,
  sprintf("lc025_class_area_unmasked_mean_%d-%d.tif", min(years), max(years))
)
retained_path <- file.path(
  out_dir,
  sprintf(
    "lc025_class_area_retained_%s_mean_%d-%d.tif",
    scenario_token,
    min(years),
    max(years)
  )
)
qa_path <- file.path(
  qa_dir,
  sprintf("lc025_class_area_weight_qa_%s.csv", scenario_token)
)

combined_mask_file <- combined_mask_path(mask_source, run_tag)
area_005_path <- here("src", "area_0p05_validdomain_km2.nc")

required_files <- c(vapply(years, input_file, character(1)), combined_mask_file, area_005_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

output_files <- c(unmasked_path, retained_path, qa_path)
outputs_current <- all(file.exists(output_files)) &&
  min(file.mtime(output_files)) >= max(file.mtime(required_files))
if (outputs_current) {
  message("Class-area weights are current; skipping calculation.")
  quit(save = "no", status = 0L)
}

area_005 <- rast(area_005_path)[[1]]
drop_mask <- rast(combined_mask_file)[[1]]
compareGeom(area_005, drop_mask, stopOnError = TRUE)
if (nrow(area_005) != 3600L || ncol(area_005) != 7200L) {
  stop("Expected a 3600 x 7200 0.05-degree area raster", call. = FALSE)
}

# terra cell values are in row-major order. Wide matrices therefore preserve
# the north-to-south row order used by the source NetCDF files.
area_values <- values(area_005, dataframe = FALSE)
drop_values <- values(drop_mask, dataframe = FALSE)
area_matrix <- matrix(
  area_values,
  nrow = nrow(area_005),
  ncol = ncol(area_005),
  byrow = TRUE
)
retained_matrix <- area_matrix
drop_matrix <- matrix(
  drop_values,
  nrow = nrow(area_005),
  ncol = ncol(area_005),
  byrow = TRUE
)
retained_matrix[is.finite(drop_matrix) & drop_matrix == 1] <- NA_real_

Rcpp::cppFunction('
NumericMatrix aggregate_class_areas(
    IntegerMatrix source,
    IntegerVector lookup,
    NumericMatrix valid_area,
    NumericMatrix retained_area,
    int class_count) {
  const int source_columns = source.nrow();
  const int source_rows = source.ncol();
  const int source_per_fine = 18;
  const int fine_per_target = 5;

  if (source_columns % source_per_fine != 0 ||
      source_rows % source_per_fine != 0) {
    stop("Source dimensions are not divisible by the 0.05-degree aggregation factor");
  }

  const int fine_columns = source_columns / source_per_fine;
  const int fine_rows = source_rows / source_per_fine;
  const int target_columns = fine_columns / fine_per_target;
  const int target_rows = fine_rows / fine_per_target;
  const int fine_cells = fine_rows * fine_columns;
  const int target_cells = target_rows * target_columns;
  if (fine_rows % fine_per_target != 0) {
    stop("Fine-grid row count is not divisible by the 0.25-degree aggregation factor");
  }
  if (class_count > 32) {
    stop("class_count exceeds fixed local buffer");
  }

  if (valid_area.nrow() != fine_rows || valid_area.ncol() != fine_columns ||
      retained_area.nrow() != fine_rows || retained_area.ncol() != fine_columns) {
    stop("Fine-grid area block has unexpected dimensions");
  }

  // Each 0.05-degree cell is independent. Calculate its class areas in
  // parallel, then aggregate the 5 x 5 fine cells to 0.25 degrees in a short
  // sequential pass. Each thread writes to a unique fine-cell column.
  NumericMatrix fine_output(2 * class_count, fine_cells);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int fine_cell = 0; fine_cell < fine_cells; ++fine_cell) {
    const int fine_row = fine_cell / fine_columns;
    const int fine_column = fine_cell % fine_columns;
    const int row_start = fine_row * source_per_fine;
    const int column_start = fine_column * source_per_fine;

    double counts[32] = {0.0};
    double total = 0.0;
    for (int row = row_start; row < row_start + source_per_fine; ++row) {
      for (int column = column_start;
           column < column_start + source_per_fine;
           ++column) {
        const int value = source(column, row);
        if (value < 0 || value >= lookup.size()) continue;
        const int class_id = lookup[value];
        if (class_id == 0) continue;
        counts[class_id - 1] += 1.0;
        total += 1.0;
      }
    }
    if (total <= 0.0) continue;

    const double area_unmasked = valid_area(fine_row, fine_column);
    const double area_retained = retained_area(fine_row, fine_column);
    const bool use_unmasked = R_finite(area_unmasked) && area_unmasked > 0.0;
    const bool use_retained = R_finite(area_retained) && area_retained > 0.0;
    if (!use_unmasked && !use_retained) continue;

    for (int class_id = 0; class_id < class_count; ++class_id) {
      const double fraction = counts[class_id] / total;
      if (use_unmasked) {
        fine_output(class_id, fine_cell) = fraction * area_unmasked;
      }
      if (use_retained) {
        fine_output(class_count + class_id, fine_cell) =
          fraction * area_retained;
      }
    }
  }

  NumericMatrix output(2 * class_count, target_cells);
  for (int fine_cell = 0; fine_cell < fine_cells; ++fine_cell) {
    const int fine_row = fine_cell / fine_columns;
    const int fine_column = fine_cell % fine_columns;
    const int target_row = fine_row / fine_per_target;
    const int target_column = fine_column / fine_per_target;
    const int target_cell = target_row * target_columns + target_column;
    for (int class_id = 0; class_id < 2 * class_count; ++class_id) {
      output(class_id, target_cell) += fine_output(class_id, fine_cell);
    }
  }
  return output;
}
', plugins = "openmp")

class_count <- length(parent_classes)
unmasked_sum <- array(0, dim = c(class_count, 720L, 1440L))
retained_sum <- array(0, dim = c(class_count, 720L, 1440L))

target_batch_size <- 5L

for (year in years) {
  message("Processing land-cover class areas for ", year)
  nc <- ncdf4::nc_open(input_file(year))
  expected_size <- c(129600L, 64800L, 1L)
  if (!identical(nc$var$lccs_class$size, expected_size)) {
    ncdf4::nc_close(nc)
    stop("Unexpected lccs_class dimensions in ", input_file(year))
  }

  for (target_start in seq.int(1L, 720L, by = target_batch_size)) {
    target_rows <- target_start:min(720L, target_start + target_batch_size - 1L)
    batch_rows <- length(target_rows)
    source_start <- (target_start - 1L) * 90L + 1L
    source <- ncdf4::ncvar_get(
      nc,
      "lccs_class",
      start = c(1L, source_start, 1L),
      count = c(129600L, 90L * batch_rows, 1L),
      collapse_degen = TRUE,
      raw_datavals = TRUE
    )
    source <- matrix(
      as.integer(source),
      nrow = 129600L,
      ncol = 90L * batch_rows
    )

    fine_start <- (target_start - 1L) * 5L + 1L
    fine_rows <- fine_start:(fine_start + 5L * batch_rows - 1L)
    row_result <- aggregate_class_areas(
      source,
      lookup,
      area_matrix[fine_rows, , drop = FALSE],
      retained_matrix[fine_rows, , drop = FALSE],
      class_count
    )

    for (batch_row in seq_len(batch_rows)) {
      output_columns <- ((batch_row - 1L) * 1440L + 1L):(batch_row * 1440L)
      target_row <- target_rows[batch_row]
      unmasked_sum[, target_row, ] <-
        unmasked_sum[, target_row, ] +
        row_result[seq_len(class_count), output_columns, drop = FALSE]
      retained_sum[, target_row, ] <-
        retained_sum[, target_row, ] +
        row_result[
          class_count + seq_len(class_count),
          output_columns,
          drop = FALSE
        ]
    }

    target_end <- max(target_rows)
    if (target_end %% 120L == 0L || target_end == 720L) {
      message("  completed target rows 1-", target_end)
    }
  }
  ncdf4::nc_close(nc)
  gc(verbose = FALSE)
}

unmasked_mean <- unmasked_sum / length(years)
retained_mean <- retained_sum / length(years)
rm(unmasked_sum, retained_sum, area_matrix, retained_matrix, drop_matrix)
gc(verbose = FALSE)

write_class_area <- function(x, path) {
  value_matrix <- vapply(
    seq_along(parent_classes),
    function(class_id) as.vector(t(x[class_id, , ])),
    numeric(720L * 1440L)
  )
  out <- rast(
    nrows = 720L,
    ncols = 1440L,
    nlyrs = length(parent_classes),
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 90,
    crs = "EPSG:4326"
  )
  values(out) <- value_matrix
  names(out) <- paste0("lc_", parent_classes)
  writeRaster(
    out,
    path,
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
  )
  out
}

unmasked_raster <- write_class_area(unmasked_mean, unmasked_path)
retained_raster <- write_class_area(retained_mean, retained_path)

unmasked_by_class <- as.numeric(global(unmasked_raster, "sum", na.rm = TRUE)[, 1])
retained_by_class <- as.numeric(global(retained_raster, "sum", na.rm = TRUE)[, 1])
qa <- tibble(
  lc_id = parent_classes,
  unmasked_class_area_km2 = unmasked_by_class,
  retained_class_area_km2 = retained_by_class,
  retained_pct = 100 * retained_by_class / unmasked_by_class
)
write_csv(
  qa,
  qa_path
)

base_area_total <- sum(area_values, na.rm = TRUE)
retained_area_total <- sum(
  ifelse(is.finite(drop_values) & drop_values == 1, NA_real_, area_values),
  na.rm = TRUE
)
message("Unmasked class-area closure: ", sum(unmasked_by_class) / base_area_total)
message("Retained class-area closure: ", sum(retained_by_class) / retained_area_total)
message("Wrote ", unmasked_path)
message("Wrote ", retained_path)
