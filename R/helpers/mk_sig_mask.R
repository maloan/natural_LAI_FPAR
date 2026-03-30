
#!/usr/bin/env Rscript
# ==============================================================================
# mk_sig_mask.R — Pixel-wise MK significance mask (p < alpha) for multi-layer raster
#
# Usage:
#   Rscript mk_sig_mask.R in_ts.nc out_mask.nc alpha [mask.nc]
#
# Notes:
#   - Reads the input as a SpatRaster with many layers (e.g., Band1_1..Band1_43).
#   - If mask.nc is provided (0/1 single layer), MK is computed only where mask==1.
#   - Output is a single-layer 0/1 raster named "mk_sig".
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
})

args <- commandArgs(trailingOnly = TRUE)
in_nc <- args[1]
out_nc <- args[2]
alpha <- as.numeric(args[3])
mask_nc <- if (length(args) >= 4) args[4] else NA_character_

r <- rast(in_nc) # multi-layer time series (Band1_1 ..)

if (!is.na(mask_nc)) {
  m <- rast(mask_nc)[[1]]
  if (!compareGeom(r[[1]], m, stopOnError = FALSE)) {
    stop("Mask geometry does not match input series: ", mask_nc)
  }
  r <- ifel(m == 1, r, NA)
}

mk_pvalue <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3) {
    return(NA_real_)
  }
  if (all(x == x[1])) {
    return(1.0)
  }

  summ <- 0L
  for (i in 1:(n - 1)) {
    d <- x[(i + 1):n] - x[i]
    summ <- summ + sum(sign(d))
  }

  tab <- table(x)
  tie_term <- sum(tab * (tab - 1) * (2 * tab + 5))
  varS <- (n * (n - 1) * (2 * n + 5) - tie_term) / 18
  if (!is.finite(varS) || varS <= 0) {
    return(1.0)
  }

  z <- if (summ > 0) (summ - 1) / sqrt(varS) else if (summ < 0) (summ + 1) / sqrt(varS) else 0
  as.numeric(stats::pnorm(-abs(z)) * 2)
}

mk_sig <- function(x) {
  p <- mk_pvalue(x)
  if (!is.finite(p)) {
    return(0)
  } # numeric scalar
  as.numeric(p < alpha) # numeric scalar 0/1
}

sig <- app(r, mk_sig)
names(sig) <- "mk_sig"

writeCDF(
  sig, out_nc,
  overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=DEFLATE"))
)

cat("Layers read:", nlyr(r), "\n")
cat("Wrote:", out_nc, "\n")
