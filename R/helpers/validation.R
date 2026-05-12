## =============================================================================
## validation.R — Input and output validation helpers
##
## Use these functions in analysis scripts to validate inputs are well-formed,
## dimensions match expectations, and outputs contain valid data.
## =============================================================================

# Validate raster dimensions match expected 0.25° global grid
validate_raster_dims_025 <- function(r, name = "raster") {
  expected_rows <- 720L
  expected_cols <- 1440L

  if (nrow(r) != expected_rows || ncol(r) != expected_cols) {
    stop(
      sprintf(
        "%s has wrong dimensions: got %d×%d, expected %d×%d (0.25° global grid)",
        name, nrow(r), ncol(r), expected_rows, expected_cols
      )
    )
  }

  invisible(r)
}

# Validate two rasters have compatible geometry
validate_geometry <- function(r1, r2, name1 = "raster1", name2 = "raster2") {
  if (!terra::compareGeom(r1, r2, stopOnError = FALSE)) {
    stop(
      sprintf(
        "Geometry mismatch between %s and %s:\n%s: %d×%d\n%s: %d×%d",
        name1, name2,
        name1, nrow(r1), ncol(r1),
        name2, nrow(r2), ncol(r2)
      )
    )
  }

  invisible(r1)
}

# Validate area weights are reasonable
validate_area_weights <- function(area, name = "area") {
  area_vals <- terra::values(area, dataframe = FALSE)

  # Check for positive finite values
  ok <- is.finite(area_vals) & area_vals > 0
  n_invalid <- sum(!ok)

  if (n_invalid > 0) {
    pct_invalid <- 100 * n_invalid / length(area_vals)
    warning(
      sprintf(
        "%s: %.1f%% pixels are invalid (NaN, Inf, or zero)",
        name, pct_invalid
      )
    )
  }

  # Check sum makes sense (should be ~150 Mkm²)
  area_total <- sum(area_vals[ok], na.rm = TRUE)
  if (area_total < 100e6 || area_total > 200e6) {
    warning(
      sprintf(
        "%s: total area is %.1f Mkm² (expected ~150 Mkm²)",
        name, area_total / 1e6
      )
    )
  }

  invisible(area)
}

# Validate trend raster values are in reasonable range
validate_trend_values <- function(r, var = "LAI", trend_type = "absolute", name = "trend") {
  vals <- terra::values(r, dataframe = FALSE)
  ok <- is.finite(vals)

  if (sum(ok) == 0) {
    stop(sprintf("%s: no finite values found", name))
  }

  vals_ok <- vals[ok]

  # Define expected ranges by variable and type
  ranges <- list(
    LAI = list(absolute = c(-0.1, 0.1), relative = c(-0.5, 0.5)),
    FPAR = list(absolute = c(-0.05, 0.05), relative = c(-0.5, 0.5))
  )

  if (!var %in% names(ranges)) {
    warning(sprintf("Unknown variable: %s, skipping range check", var))
    return(invisible(r))
  }

  expected_range <- ranges[[var]][[trend_type]]

  # Check percentiles
  p01 <- quantile(vals_ok, 0.01, na.rm = TRUE, type = 7)
  p99 <- quantile(vals_ok, 0.99, na.rm = TRUE, type = 7)

  if (p01 < expected_range[1] * 1.5 || p99 > expected_range[2] * 1.5) {
    warning(
      sprintf(
        "%s (%s %s): 1st-99th percentile range [%.4f, %.4f], expected ~[%.4f, %.4f]",
        name, var, trend_type, p01, p99, expected_range[1], expected_range[2]
      )
    )
  }

  invisible(r)
}

# Validate output CSV has required columns
validate_output_csv <- function(df, required_cols = NULL, name = "dataframe") {
  if (is.null(required_cols)) {
    # Default columns for trend summaries
    required_cols <- c("scenario", "mean", "ci_lower", "ci_upper", "n_pixels", "area_km2")
  }

  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "%s missing required columns: %s",
        name, paste(missing_cols, collapse = ", ")
      )
    )
  }

  invisible(df)
}

# Validate CSV metadata (units, methods documented)
validate_csv_metadata <- function(file_path) {
  lines <- readLines(file_path, n = 10)
  has_metadata <- any(grepl("^#", lines))

  if (!has_metadata) {
    warning(sprintf("%s: no metadata comments found (add with write_csv(..., comment='...'))", file_path))
  }

  invisible(file_path)
}

# Validate all files in an output directory
validate_output_dir <- function(dir_path, expected_files = NULL) {
  if (!dir.exists(dir_path)) {
    stop(sprintf("Output directory does not exist: %s", dir_path))
  }

  files <- list.files(dir_path)

  if (length(files) == 0) {
    stop(sprintf("Output directory is empty: %s", dir_path))
  }

  if (!is.null(expected_files)) {
    missing <- setdiff(expected_files, files)
    if (length(missing) > 0) {
      warning(
        sprintf(
          "Missing expected files in %s: %s",
          dir_path, paste(missing, collapse = ", ")
        )
      )
    }
  }

  invisible(dir_path)
}
