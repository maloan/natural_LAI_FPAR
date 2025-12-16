## =============================================================================
## 22_pixel_trend_significance.R
## Pixel-wise trend significance (OLS + FDR), 0.25°
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(here)
  library(readr)
})

ROOT   <- here::here()
DIR025 <- file.path(ROOT, "analysis", "unmasked", "0p25")
OUTDIR <- file.path(ROOT, "analysis", "pixel_significance")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

VARS    <- c("LAI", "FPAR")
METRICS <- c("yearmean", "yearmax")

ALPHA_P <- 0.05
ALPHA_Q <- 0.10

pixel_lm_stats <- function(v, years) {
  ok <- is.finite(v)
  if (sum(ok) < 10) {
    return(c(slope = NA_real_, p = NA_real_))
  }
  fit <- lm(v[ok] ~ years[ok])
  c(
    slope = coef(fit)[2],
    p     = summary(fit)$coefficients[2, 4]
  )
}



for (var in VARS) {
  for (met in METRICS) {

    message("Processing: ", var, " / ", met)

    f_nc <- file.path(DIR025,
                      sprintf("%s_georef_%s_0p25.nc", var, met))
    if (!file.exists(f_nc)) next

    r <- rast(f_nc)

    # ------------------------------------------------------------------
    # Per-pixel OLS: slope + nominal p-value
    # ------------------------------------------------------------------
    yrs <- 1982:(1982 + nlyr(r) - 1L)
    r_stats <- terra::app(r, \(x) pixel_lm_stats(x, yrs))
    names(r_stats) <- c("slope_per_year", "p_value")

    # ------------------------------------------------------------------
    # Save slope and p-value rasters
    # ------------------------------------------------------------------
    writeCDF(
      r_stats,
      file.path(
        OUTDIR,
        sprintf("%s_%s_pixel_slope_p.nc", var, met)
      ),
      overwrite = TRUE
    )

    # ------------------------------------------------------------------
    # FDR correction (Benjamini–Hochberg)
    # ------------------------------------------------------------------
    df <- as.data.frame(r_stats, xy = FALSE, na.rm = FALSE)

    valid <- is.finite(df$p_value)
    qvals <- rep(NA_real_, nrow(df))
    qvals[valid] <- p.adjust(df$p_value[valid], method = "BH")

    # back to raster
    r_q <- setValues(r_stats[[1]], qvals)
    names(r_q) <- "q_value"

    writeCDF(
      r_q,
      file.path(
        OUTDIR,
        sprintf("%s_%s_pixel_q.nc", var, met)
      ),
      overwrite = TRUE
    )

    r_sig_p <- r_stats[["p_value"]] < ALPHA_P
    r_sig_q <- r_q < ALPHA_Q

    names(r_sig_p) <- "sig_p_0p05"
    names(r_sig_q) <- "sig_q_0p10"

    writeCDF(
      c(r_sig_p, r_sig_q),
      file.path(
        OUTDIR,
        sprintf("%s_%s_pixel_significance_masks.nc", var, met)
      ),
      overwrite = TRUE
    )

  }
}

message("Finished pixel-wise trend significance.")

