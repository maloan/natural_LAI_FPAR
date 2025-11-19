## =============================================================================
# names.R — Centralized naming for derived raster products
#
# Purpose
#   Consistent, validated filename templates for categorical and fractional
#   land-cover rasters. Ensures stable conventions for dataset tokens,
#   component names, and resolution identifiers.
#
# Notes
#   - Base R only.
#   - Enforces valid years (YYYY).
#   - Normalises arbitrary dataset/kind strings.
#   - Resolution fixed to "0p05".
## =============================================================================


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

.assert_year <- function(year) {
  if (!is.numeric(year) || length(year) != 1L || is.na(year))
    stop("`year` must be a single numeric value.", call. = FALSE)

  y <- as.integer(year)
  if (y < 1800L || y > 2200L)
    stop("`year` invalid: ", y, " (expected 4-digit year).", call. = FALSE)

  y
}

.norm_token <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x))
    stop("Token must be a non-empty scalar string.", call. = FALSE)

  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_", "", sub("_$", "", x))

  if (!nzchar(x))
    stop("Token normalised to empty string.", call. = FALSE)

  x
}

.res_tok <- function(dx = 0.05) {
  "0p05"
}


# ------------------------------------------------------------------------------
# Input Data
# ------------------------------------------------------------------------------

# Categorical yearstack (0.05°)
name_stack_cat <- function(dataset = c("cci", "landsat")) {
  ds <- match.arg(dataset)
  sprintf("%s_cat_yearstack_%s.tif", .norm_token(ds), .res_tok())
}

# Fractional cover component (CCI-derived), per year
name_frac <- function(kind, year) {
  k <- .norm_token(kind)
  y <- .assert_year(year)
  sprintf("cci_frac_%s_%d_%s.tif", k, y, .res_tok())
}

# Single-year categorical file (any dataset)
name_cat_year <- function(dataset, year) {
  ds <- .norm_token(dataset)
  y  <- .assert_year(year)
  sprintf("%s_cat_%d_%s.tif", ds, y, .res_tok())
}
