## =============================================================================
# names.R — Centralized naming conventions for derived raster products
#
# Purpose
#   Provide validated, consistent filename templates for categorical and
#   fractional land-cover rasters used across the workflow.
#
# Notes
#   - Base R only.
#   - Enforces YYYY, normalizes dataset/kind tokens, and fixes resolution token.
## =============================================================================

# ---- internals ---------------------------------------------------------------

.assert_year <- function(year) {
  if (!is.numeric(year) || length(year) != 1L || is.na(year)) {
    stop("`year` must be a single numeric value.", call. = FALSE)
  }
  y <- as.integer(year)
  if (y < 1800L || y > 2200L) {
    stop("`year` looks invalid: ", y, " (expected 4-digit year).", call. = FALSE)
  }
  y
}

.norm_token <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("Token must be a non-empty scalar string.", call. = FALSE)
  }
  # keep [a-z0-9_]; convert others to underscores, trim doubles/edges
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_", "", sub("_$", "", x))
  if (!nzchar(x)) {
    stop("Token normalised to empty string.", call. = FALSE)
  }
  x
}

.res_tok <- function(dx = 0.05) {
  # Use fixed convention "0p05" for 0.05°
  # If you later generalize, format here.
  "0p05"
}

# ---- public API --------------------------------------------------------------

# Yearstack filename for categorical land cover (0.05°)
# dataset: one of "cci", "landsat" (extensible; validated as a token)
name_stack_cat <- function(dataset = c("cci", "landsat")) {
  ds <- match.arg(dataset)
  sprintf("%s_cat_yearstack_%s.tif", .norm_token(ds), .res_tok())
}

# Fractional cover filename by component and year (CCI-derived fractions)
# kind: e.g., "grass", "bare", "urban", "cropland", "frac_fused"
name_frac <- function(kind, year) {
  k <- .norm_token(kind)
  y <- .assert_year(year)
  sprintf("cci_frac_%s_%d_%s.tif", k, y, .res_tok())
}

# Single-year categorical raster filename by dataset
# dataset: e.g., "cci", "landsat", "glc_fcs30d" (tokenized)
name_cat_year <- function(dataset, year) {
  ds <- .norm_token(dataset)
  y  <- .assert_year(year)
  sprintf("%s_cat_%d_%s.tif", ds, y, .res_tok())
}
