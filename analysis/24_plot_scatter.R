#!/usr/bin/env Rscript
## =============================================================================
# 24_plot_scatter.R — Scatterplots for slope comparisons
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(here)
})

# ---------------------------------------------------------------------
# Root and output paths
# ---------------------------------------------------------------------
ROOT <- here()
OUTP <- here("analysis", "trend_plots", "scatter")
dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# discover input slope files
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

slopes <- lapply(files, function(f){
  list(info = info(f), file = f, r = rast(f))
})

# ---------------------------------------------------------------------
# prepare dataframe for scatter
# ---------------------------------------------------------------------
df <- function(r1, r2){
  xy <- cbind(values(r1), values(r2))
  d  <- as.data.frame(xy)
  d  <- d[complete.cases(d), ]
  colnames(d) <- c("x","y")
  d
}

# ---------------------------------------------------------------------
# scatter plotting helper
# ---------------------------------------------------------------------
plot_sc <- function(d, ttl, of){
  g <- ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.2, size = 0.4, color = "#1f78b4") +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    coord_fixed() +
    labs(title = ttl, x = "Variant A", y = "Variant B") +
    theme_minimal(base_size = 14)
  ggsave(of, g, width = 7, height = 6)
}

# ---------------------------------------------------------------------
# masked vs unmasked scatterplots
# ---------------------------------------------------------------------
unmasked <- slopes[sapply(slopes, function(x) x$info$mask == "UNMASKED")]
masked   <- slopes[sapply(slopes, function(x) x$info$mask != "UNMASKED")]

for (um in unmasked){
  for (mk in masked){
    if (um$info$var != mk$info$var) next

    d <- df(um$r, mk$r)
    ttl <- sprintf("%s: masked (%s τ=%s) vs unmasked",
                   mk$info$var, mk$info$mask, mk$info$tau)

    ofile <- file.path(
      OUTP,
      sprintf("scatter_%s_masked_%s_%s_vs_unmasked.png",
              mk$info$var, mk$info$mask, mk$info$tau)
    )

    plot_sc(d, ttl, ofile)
  }
}

cat("Scatterplots saved in:", OUTP, "\n")
