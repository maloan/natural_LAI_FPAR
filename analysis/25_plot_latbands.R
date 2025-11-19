#!/usr/bin/env Rscript
## =============================================================================
# 25_plot_latbands.R — Latitudinal mean slope profiles
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(here)
})

# ---------------------------------------------------------------------
# Root + output directory
# ---------------------------------------------------------------------
ROOT <- here()
OUTP <- here("analysis", "trend_plots", "latbands")
dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# discover slope files
# ---------------------------------------------------------------------
files <- list.files(
  ROOT,
  pattern = "trend_slope_yearmean_decade.*\\.nc$",
  full.names = TRUE,
  recursive = TRUE
)

info <- function(f){
  var  <- if (grepl("FPAR", f)) "FPAR" else "LAI"
  mask <- if (grepl("CCI", f)) "CCI" else if (grepl("GLC", f)) "GLC" else "UNMASKED"
  tau  <- regmatches(f, regexpr("tau_0\\.[0-9]+", f))
  if (length(tau) == 0) tau <- "none"
  list(var = var, mask = mask, tau = tau)
}

# ---------------------------------------------------------------------
# compute latitudinal mean trend
# ---------------------------------------------------------------------
lat_summary <- function(r){
  xy   <- crds(r, df = TRUE)
  vals <- values(r)
  d    <- cbind(xy, slope = vals)
  d    <- d[complete.cases(d), ]
  d$latbin <- cut(d$y, breaks = seq(-90, 90, 5), include.lowest = TRUE)
  d %>% group_by(latbin) %>% summarise(mean_slope = mean(slope))
}

# ---------------------------------------------------------------------
# loop through all files
# ---------------------------------------------------------------------
for (f in files){
  inf <- info(f)
  r   <- rast(f)
  df  <- lat_summary(r)

  outdir <- file.path(OUTP, inf$var, inf$mask)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  of <- file.path(
    outdir,
    sprintf("latband_%s_%s_%s.png", inf$var, inf$mask, inf$tau)
  )

  g <- ggplot(df, aes(x = latbin, y = mean_slope)) +
    geom_col(fill = hcl.colors(5, "batlow")[3]) +
    coord_flip() +
    labs(
      title = sprintf("%s %s %s", inf$var, inf$mask, inf$tau),
      x = "Latitude band",
      y = "Slope per decade"
    ) +
    theme_minimal(base_size = 14)

  ggsave(of, g, width = 8, height = 7)
}

cat("Latitude profiles saved to:", OUTP, "\n")
