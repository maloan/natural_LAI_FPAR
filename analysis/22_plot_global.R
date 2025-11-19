#!/usr/bin/env Rscript
## =============================================================================
# 22_plot_global.R — Global trend maps (slope per decade, 0.25°)
# Palette: batlow (scientific CDT)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(here)
})

# --- root and output -----------------------------------------------------------
ROOT <- here()

OUTP <- here("analysis", "trend_plots", "global")
dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)

# --- coastline -----------------------------------------------------------------
coast <- NULL
coast_path <- here("src", "ne_110m_coastline.gpkg")
if (file.exists(coast_path)) coast <- read_sf(coast_path)

# --- palettes ------------------------------------------------------------------
pal      <- hcl.colors(100, "batlow")
pal_div  <- hcl.colors(100, "Blue-Red 3")

# --- plotting helper -----------------------------------------------------------
plot_global <- function(r, title, file, lim = c(-0.05, 0.05)) {
  png(file, width = 1800, height = 900, res = 150)
  par(mai = c(0.6, 0.8, 0.7, 0.8))
  plot(r, col = pal_div, range = lim, axes = FALSE, main = title)
  if (!is.null(coast)) plot(coast$geom, add = TRUE, lwd = 0.4)
  dev.off()
}

# --- find input trend files ----------------------------------------------------
files <- list.files(
  ROOT,
  pattern = "trend_slope_yearmean_decade.*\\.nc$",
  full.names = TRUE,
  recursive = TRUE
)

# --- main loop -----------------------------------------------------------------
for (f in files) {

  r <- rast(f)

  var  <- if (grepl("FPAR", f)) "FPAR" else "LAI"
  mask <- if (grepl("CCI", f)) "CCI" else if (grepl("GLC", f)) "GLC" else "UNMASKED"

  tau <- regmatches(f, regexpr("tau_0\\.[0-9]+", f))
  if (length(tau) == 0) tau <- "none"

  outdir <- file.path(OUTP, var, mask, tau)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  ofile <- file.path(
    outdir,
    sprintf("%s_%s_%s_global_trend.png", var, mask, tau)
  )

  title <- sprintf(
    "%s — %s — τ=%s\nTrend (yearmean per decade)",
    var, mask, tau
  )

  plot_global(r, title, ofile)
}

cat("Global trend plots saved to:", OUTP, "\n")
