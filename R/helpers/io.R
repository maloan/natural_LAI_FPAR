
## =============================================================================
## io.R — Input/output utilities for reproducible GeoTIFF and metadata handling
## =============================================================================
# ------------------------------------------------------------------------------
# GDAL creation options
# ------------------------------------------------------------------------------

gdal_co <- function(predictor = 2L, speed_over_size = FALSE) {
  base <- c("TILED=YES", "BIGTIFF=IF_SAFER")
  comp <- if (isTRUE(speed_over_size)) {
    "COMPRESS=LZW"
  } else {
    "COMPRESS=DEFLATE"
  }
  c(base, comp, sprintf("PREDICTOR=%d", predictor), "NUM_THREADS=ALL_CPUS")
}

# Legacy wrappers for backward compatibility
gdal_co_int <- function(speed_over_size = FALSE) {
  gdal_co(predictor = 2L, speed_over_size = speed_over_size)
}

gdal_co_f32 <- function(speed_over_size = FALSE) {
  gdal_co(predictor = 3L, speed_over_size = speed_over_size)
}

# ------------------------------------------------------------------------------
# Write-option
# ------------------------------------------------------------------------------
wopt_byte <- function(speed_over_size = FALSE, na = 255L) {
  list(
    datatype = "INT1U",
    NAflag   = as.integer(na),
    gdal     = gdal_co_int(speed_over_size)
  )
}

wopt_int <- function(speed_over_size = FALSE, na = NA_integer_) {
  out <- list(
    datatype = "INT2U",
    gdal     = gdal_co_int(speed_over_size)
  )
  if (!is.na(na)) out$NAflag <- as.integer(na)
  out
}

wopt_f32 <- function(speed_over_size = FALSE, na = -9999) {
  list(
    datatype = "FLT4S",
    NAflag   = as.numeric(na),
    gdal     = gdal_co_f32(speed_over_size)
  )
}
