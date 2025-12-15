
## =============================================================================
## 20_delta_yearmax_yearmean.R — Spatial Δ-trends (yearmax – yearmean) at 0.25°
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(sf)
  library(rnaturalearth)
  library(scico)
})

ROOT     <- here::here()

# Unmasked trends (from analysis/unmasked/0p25)
DIR_UNM  <- file.path(ROOT, "analysis", "unmasked", "0p25")

OUTDIR   <- file.path(ROOT, "analysis", "delta_trends")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

VARS <- c("LAI", "FPAR")

coast <- ne_coastline(scale = 110, returnclass = "sf")

for (var in VARS) {
  ## --------------------------------------------------------------------------
  ## Load yearmax and yearmean trend maps
  ## --------------------------------------------------------------------------
  f_yearmax <- file.path(DIR_UNM,
                         sprintf("%s_georef_yearmax_trend_slope_peryear_0p25.nc", var))

  f_yearmean <- file.path(DIR_UNM,
                          sprintf("%s_georef_yearmean_trend_slope_peryear_0p25.nc", var))

  if (!file.exists(f_yearmax) | !file.exists(f_yearmean)) {
    warning("Missing trend inputs for: ", var)
    next
  }

  r_yearmax  <- rast(f_yearmax)
  r_yearmean <- rast(f_yearmean)

  ## --------------------------------------------------------------------------
  ## Δ = yearmax – yearmean
  ## --------------------------------------------------------------------------
  r_delta <- r_yearmax - r_yearmean

  ## --------------------------------------------------------------------------
  ## Save NetCDF
  ## --------------------------------------------------------------------------
  out_nc <- file.path(OUTDIR,
                      sprintf("DELTA_%s_yearmax_minus_yearmean_0p25.nc", var))
  writeCDF(r_delta, out_nc, overwrite = TRUE)

  ## --------------------------------------------------------------------------
  ## Quicklook PNG
  ## --------------------------------------------------------------------------
  df <- as.data.frame(r_delta, xy = TRUE, na.rm = FALSE)
  names(df)[1:2] <- c("x", "y")
  names(df)[3]   <- "delta"

  ## symmetric ±2σ
  lim <- 2 * sd(df$delta, na.rm = TRUE)

  p <- ggplot(df, aes(x = x, y = y, fill = delta)) +
    geom_raster() +
    geom_sf(
      data = coast,
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.2
    ) +
    scale_fill_scico(
      palette = "bam",
      name = expression(Delta == yearmax - yearmean),
      limits = c(-lim, lim)
    ) +
    coord_sf(expand = FALSE) +
    labs(title = sprintf("Δ-trend(yearmax – yearmean) — %s", var),
         subtitle = "Positive = stronger trend in seasonal peak (yearmax)") +
    theme_bw(base_size = 11)

  out_png <- file.path(OUTDIR,
                       sprintf("DELTA_%s_yearmax_minus_yearmean_0p25.png", var))

  ggsave(out_png,
         p,
         width = 8,
         height = 4.5,
         dpi = 330)
  message("Saved Δ map: ", out_png)
}

message("Finished Δ-trend computation.")

