## =============================================================================
## io.R — Input/output utilities for reproducible GeoTIFF and metadata handling
## =============================================================================
# ------------------------------------------------------------------------------
# GDAL creation options
# ------------------------------------------------------------------------------

gdal_co_int <- function(speed_over_size = FALSE) {
  base <- c("TILED=YES", "BIGTIFF=IF_SAFER")
  comp <- if (isTRUE(speed_over_size)) {
    "COMPRESS=LZW"
  } else {
    "COMPRESS=DEFLATE"
  }
  c(base, comp, "PREDICTOR=2", "NUM_THREADS=ALL_CPUS")
}

gdal_co_f32 <- function(speed_over_size = FALSE) {
  base <- c("TILED=YES", "BIGTIFF=IF_SAFER")
  comp <- if (isTRUE(speed_over_size)) {
    "COMPRESS=LZW"
  } else {
    "COMPRESS=DEFLATE"
  }
  c(base, comp, "PREDICTOR=3", "NUM_THREADS=ALL_CPUS")
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
