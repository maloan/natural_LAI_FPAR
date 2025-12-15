## =============================================================================
## 21_trend_significance_maps.R — Pixel-wise trend significance (0.25° annual)
##   - unmasked + masked (CCI/GLC × τ)
##   - AR(1) effective sample size corrected p-values
##   - optional BH-FDR q-values
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(scico)
  library(sf)
  library(rnaturalearth)
  library(scales)
})

ROOT <- here::here()

# ------------------------------------------------------------------------------
# USER CONFIG
# ------------------------------------------------------------------------------
VARS    <- c("LAI", "FPAR")
METRICS <- c("yearmean", "yearmax", "yearamp")

MASKS <- c("CCI", "GLC")
TAU   <- c("tau_0.05", "tau_0.1", "tau_0.2")

IN_UNMASK <- file.path(ROOT, "analysis", "unmasked", "0p25")
OUTDIR    <- file.path(ROOT, "analysis", "trend_significance")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ALPHA      <- 0.05      # significance level
USE_FDR    <- TRUE      # also compute BH q-values and mask on q<ALPHA
SD_K       <- 2         # map clamping
MIN_N      <- 12        # minimum non-NA years to fit trend

NCORES <- as.integer(Sys.getenv("NCORES", "1"))
terraOptions(progress = 1)

coast   <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
lab_deg <- scales::label_number(suffix = "°")

# ------------------------------------------------------------------------------
# Per-pixel trend p-value (AR1-effective sample size correction)
# ------------------------------------------------------------------------------
trend_p_ar1 <- function(y) {
  # y: numeric vector (time series for one pixel)
  ok <- is.finite(y)
  n  <- sum(ok)
  if (n < MIN_N) return(c(NA_real_, NA_real_, NA_real_))  # slope, p, neff

  yy <- y[ok]
  tt <- seq_len(n)

  fit <- try(lm(yy ~ tt), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(NA_real_, NA_real_, NA_real_))

  slope <- unname(coef(fit)[2])

  res <- residuals(fit)
  if (length(res) < 3) return(c(slope, NA_real_, NA_real_))

  # lag-1 autocorrelation of residuals
  r1 <- suppressWarnings(stats::cor(res[-length(res)], res[-1], use = "complete.obs"))
  if (!is.finite(r1)) r1 <- 0
  r1 <- max(min(r1, 0.99), -0.99)

  # effective sample size (approx)
  neff <- n * (1 - r1) / (1 + r1)
  neff <- max(neff, 3.1)  # keep df > 1

  # inflate SE by sqrt(n / neff)
  se_ols <- summary(fit)$coefficients[2, 2]
  se_adj <- se_ols * sqrt(n / neff)

  tval <- slope / se_adj
  df   <- neff - 2
  pval <- 2 * pt(-abs(tval), df = df)

  c(slope, pval, neff)
}

# ------------------------------------------------------------------------------
# Plot helper: slope map with significance fading
# ------------------------------------------------------------------------------
plot_slope_sig <- function(df, var, metric, domain_lab, SD_K = 2) {
  sdev  <- sd(df$slope, na.rm = TRUE)
  clamp <- SD_K * sdev
  df$slope_c <- pmax(pmin(df$slope, clamp), -clamp)

  ggplot(df, aes(lon, lat)) +
    geom_raster(aes(fill = slope_c, alpha = sig)) +
    geom_sf(data = coast, inherit.aes = FALSE,
            colour = "black", linewidth = 0.2) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 60), labels = lab_deg) +
    scale_y_continuous(breaks = seq(-90, 90, 30), labels = lab_deg) +
    scale_fill_scico(
      palette = "bam",
      name    = paste0(var, " trend"),
      limits  = c(-clamp, clamp),
      oob     = scales::squish
    ) +
    scale_alpha_manual(values = c(`0` = 0.20, `1` = 1.00), guide = "none") +
    labs(
      title    = sprintf("%s %s trend (significance faded)", var, metric),
      subtitle = domain_lab,
      x = "Longitude", y = "Latitude"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank())
}

# ------------------------------------------------------------------------------
# Run one domain (unmasked or masked)
# ------------------------------------------------------------------------------
run_domain <- function(var, metric, domain, tau = NA, mask = NA) {

  if (domain == "unmasked") {
    f_data  <- file.path(IN_UNMASK, sprintf("%s_georef_%s_0p25.nc", var, metric))
    f_slope <- file.path(IN_UNMASK, sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, metric))
    domain_lab <- "Unmasked (all land)"
    out_base <- file.path(OUTDIR, "unmasked", var, metric)
  } else {
    base <- file.path(ROOT, "output", tau, "eval", sprintf("trend_%s_%s", var, mask))
    f_data  <- file.path(base, sprintf("%s_%s_0p25.nc", var, metric))
    f_slope <- file.path(base, sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, metric))
    domain_lab <- sprintf("Masked (%s, τ=%s)", mask, sub("tau_", "", tau))
    out_base <- file.path(OUTDIR, "masked", tau, paste0(var, "_", mask), metric)
  }

  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(f_data) || !file.exists(f_slope)) {
    warning("Missing inputs: ", f_data, " or ", f_slope)
    return(invisible(NULL))
  }

  r_data  <- rast(f_data)
  r_slope <- rast(f_slope)

  # compute slope + p + neff from the actual time series (not from stored slope)
  # (keeps p-values coherent even if slope file was computed differently)
  message("Computing p-values: ", domain_lab, " / ", var, " / ", metric)

  r_stats <- terra::app(r_data, fun = trend_p_ar1, cores = NCORES)
  names(r_stats) <- c("slope_fit", "pval", "neff")

  # keep your original slope raster (for consistency with earlier figures),
  # but use p-values from the fit
  r_pval <- r_stats[["pval"]]

  # significance mask
  r_sig <- ifel(r_pval < ALPHA, 1, 0)

  # optional FDR (BH) q-values
  r_qval <- NULL
  if (USE_FDR) {
    p <- values(r_pval, mat = FALSE)
    ok <- is.finite(p)
    q <- rep(NA_real_, length(p))
    q[ok] <- p.adjust(p[ok], method = "BH")
    r_qval <- r_pval
    values(r_qval) <- q
    r_sig <- ifel(r_qval < ALPHA, 1, 0)
  }

  # save rasters
  writeCDF(r_pval, file.path(out_base, sprintf("%s_%s_trend_pval_ar1_0p25.nc", var, metric)),
           overwrite = TRUE)
  writeCDF(r_sig,  file.path(out_base, sprintf("%s_%s_trend_sigmask_0p25.nc", var, metric)),
           overwrite = TRUE)
  if (!is.null(r_qval)) {
    writeCDF(r_qval, file.path(out_base, sprintf("%s_%s_trend_qval_bh_0p25.nc", var, metric)),
             overwrite = TRUE)
  }

  # map dataframe (slope from your stored slope file; sig from mask)
  df_slope <- as.data.frame(r_slope, xy = TRUE, na.rm = FALSE)
  colnames(df_slope)[1:2] <- c("lon", "lat")
  colnames(df_slope)[3]   <- "slope"

  df_sig <- as.data.frame(r_sig, xy = TRUE, na.rm = FALSE)
  colnames(df_sig)[1:2] <- c("lon", "lat")
  colnames(df_sig)[3]   <- "sig"

  df <- df_slope |>
    inner_join(df_sig, by = c("lon", "lat")) |>
    filter(is.finite(slope), is.finite(sig))

  p_map <- plot_slope_sig(df, var, metric, domain_lab, SD_K = SD_K)

  ggsave(
    file.path(out_base, sprintf("%s_%s_trend_map_sigfaded.png", var, metric)),
    p_map, width = 7.2, height = 3.8, dpi = 330
  )

  # global area fraction significant (area-weighted, by latitude)
  df_frac <- df |>
    mutate(w = cos(lat * pi / 180)) |>
    summarise(
      frac_sig = sum(w * (sig == 1), na.rm = TRUE) / sum(w, na.rm = TRUE),
      n_pix    = n()
    )

  write_csv(df_frac, file.path(out_base, "summary_significance_fraction.csv"))

  invisible(NULL)
}

# ------------------------------------------------------------------------------
# Main loops
# ------------------------------------------------------------------------------
# Unmasked
for (var in VARS) {
  for (met in METRICS) {
    run_domain(var, met, domain = "unmasked")
  }
}

# Masked
for (tau in TAU) {
  for (var in VARS) {
    for (mask in MASKS) {
      for (met in METRICS) {
        run_domain(var, met, domain = "masked", tau = tau, mask = mask)
      }
    }
  }
}

message("Finished pixel-wise trend significance.")
