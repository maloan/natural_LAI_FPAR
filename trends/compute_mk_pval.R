# ==============================================================================
# compute_mk_pval.R
# Pixelwise Mann–Kendall p-values for annual-mean (or other) time series stacks.
# Compute p-value of monotonic trend (Mann–Kendall) for each pixel in a
# SpatRaster time stack, writing a single-layer NetCDF.
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
  library(trend)
})

# ---- config
mode <- "masked"
var <- "LAI"
metric <- "yearmean"

tau <- "tau_0.1"
mask <- "CCI"

fill_value <- -9999

# ---- paths
if (mode == "unmasked") {
  f_ts <- here(
    "analysis", "unmasked", "0p25",
    sprintf("%s_georef_%s_0p25.nc", var, metric)
  )
  out_p <- here(
    "analysis", "unmasked", "0p25",
    sprintf("%s_georef_%s_trend_mk_pval_0p25.nc", var, metric)
  )
} else if (mode == "masked") {
  f_ts <- here(
    "output", tau, "eval",
    sprintf("trend_%s_%s", var, mask),
    sprintf("%s_%s_0p25.nc", var, metric)
  )
  out_p <- here(
    "output", tau, "eval",
    sprintf("trend_%s_%s", var, mask),
    sprintf("%s_%s_trend_mk_pval_0p25.nc", var, metric)
  )
} else {
  stop("mode must be 'masked' or 'unmasked'")
}

if (!file.exists(f_ts)) stop("Input time stack not found: ", f_ts)
dir.create(dirname(out_p), recursive = TRUE, showWarnings = FALSE)

message("Input : ", f_ts)
message("Output: ", out_p)

# ---- load & prep -------------------------------------------------------------
r <- terra::rast(f_ts)
r[r == fill_value] <- NA

# ---- pixel function (Mann–Kendall p-value)
mk_p_fun <- function(y) {
  ok <- is.finite(y)
  yy <- y[ok]
  n <- length(yy)
  if (n < 3) {
    return(NA_real_)
  }

  # trend::mk.test returns a list; p.value is the two-sided p-value
  out <- try(trend::mk.test(yy)$p.value, silent = TRUE)
  if (inherits(out, "try-error") || !is.finite(out)) {
    return(NA_real_)
  }
  out
}

# ---- compute & write
pval <- terra::app(r, mk_p_fun)
names(pval) <- "pval_mk"

terra::writeCDF(pval, out_p, overwrite = TRUE)
