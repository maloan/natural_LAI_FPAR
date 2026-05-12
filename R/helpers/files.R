## =============================================================================
## files.R — Filename and file-discovery helpers
## =============================================================================

#' Format decimal as zero-padded string with 'p' for decimal point
#'
#' Converts numeric value to string formatted as "0.XX" with decimal point
#' replaced by 'p' (e.g., 0.05 → "0p05"). Useful for filename generation.
#'
#' @param x Numeric value to format
#' @return Character string with 'p' replacing decimal point
#'
#' @examples
#' tok(0.05) # "0p05"
#' tok(0.1) # "0p10"
#'
tok <- function(x) gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))

#' Extract YYYY-MM from filename
#'
#' Searches filename for ISO 8601 YYYYMM pattern and returns it.
#' Falls back to basename if no match found.
#'
#' @param p File path
#' @return Character string with YYYYMM date or basename
#'
extract_ym_from_filename <- function(p) {
  s <- tools::file_path_sans_ext(basename(p))
  m <- regexpr("(19|20)\\d{2}(0[1-9]|1[0-2])", s, perl = TRUE)
  if (m[1] > 0) {
    substr(s, m[1], m[1] + attr(m, "match.length") - 1)
  } else {
    s
  }
}

#' Find exactly one file matching a pattern
#'
#' Searches directory for files matching a regex pattern and returns the path.
#' Errors if zero or multiple files match (must be exactly one).
#'
#' @param dir Directory to search
#' @param pattern Regex pattern to match filenames
#' @param label Optional description for error messages (e.g., "LAI mask")
#'
#' @return Full path to the single matching file
#'
find_one <- function(dir, pattern, label = NULL) {
  cand <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(cand) != 1L) {
    msg <- if (!is.null(label)) {
      paste0(label, ": expected 1 file in ", dir, ", found ", length(cand))
    } else {
      paste0("Expected 1 file matching '", pattern, "' in ", dir, ", found ", length(cand))
    }
    stop(msg, call. = FALSE)
  }
  cand
}

#' Filter indices within a year range
#'
#' Returns indices where corresponding values fall within [year_min, year_max].
#' Useful for subsetting year vectors or layer stacks by year window.
#'
#' @param values Numeric vector of year values
#' @param year_min Minimum year (inclusive)
#' @param year_max Maximum year (inclusive)
#'
#' @return Integer vector of indices where year_min <= values <= year_max
#'
filter_by_year_range <- function(values, year_min, year_max) {
  which(values >= year_min & values <= year_max)
}

#' Align raster to template geometry
#'
#' Checks if raster geometry matches template. If not, resamples to match.
#' Uses terra::compareGeom and terra::resample internally.
#'
#' @param r SpatRaster to align
#' @param template SpatRaster template for target geometry
#' @param method Resampling method: "bilinear", "near", etc. (terra methods)
#'
#' @return Raster aligned to template geometry
#'
align_to_template <- function(r, template, method = "bilinear") {
  if (!terra::compareGeom(r, template, stopOnError = FALSE)) {
    r <- terra::resample(r, template, method = method)
  }
  r
}

#' Stop with standardized error message
#'
#' Convenience wrapper for stop() that automatically suppresses the call stack.
#' Results in cleaner error messages for end users.
#'
#' @param ... Character strings and expressions to concatenate into error message
#'
stop_msg <- function(...) {
  stop(paste0(...), call. = FALSE)
}
