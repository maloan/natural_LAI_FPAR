#!/usr/bin/env Rscript
## =============================================================================
# 23_plot_difference.R — Difference maps for trend diagnostics
#   - masked − unmasked
#   - CCI τ comparisons
#   - GLC − CCI
# Palette: Blue-Red 3 (diverging)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(here)
})

# ---------------------------------------------------------------------
# Root + output using here()
# ---------------------------------------------------------------------
ROOT <- here()
OUTP <- here("analysis", "trend_plots", "differences")
dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)

pal_div <- hcl.colors(100, "Blue-Red 3")

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
  if (length(tau)==0) tau <- "none"
  list(var = var, mask = mask, tau = tau)
}

slopes <- lapply(files, function(f){
  list(info = info(f), file = f, r = rast(f))
})

# ---------------------------------------------------------------------
# plotting helper
# ---------------------------------------------------------------------
plot_diff <- function(r, title, ofile){
  png(ofile, width = 1800, height = 900, res = 150)
  plot(r, col = pal_div, range = c(-0.05, 0.05), axes = FALSE, main = title)
  dev.off()
}

# ---------------------------------------------------------------------
# masked − unmasked
# ---------------------------------------------------------------------
unmasked <- slopes[sapply(slopes, function(x) x$info$mask == "UNMASKED")]
masked   <- slopes[sapply(slopes, function(x) x$info$mask != "UNMASKED")]

for (um in unmasked){
  var <- um$info$var

  for (mk in masked){
    if (mk$info$var != var) next

    d  <- mk$r - um$r
    of <- file.path(
      OUTP,
      sprintf("%s_masked_%s_%s_minus_unmasked.png",
              var, mk$info$mask, mk$info$tau)
    )

    ttl <- sprintf("%s: Masked (%s τ=%s) − Unmasked",
                   var, mk$info$mask, mk$info$tau)

    plot_diff(d, ttl, of)
  }
}

# ---------------------------------------------------------------------
# CCI τ comparisons
# ---------------------------------------------------------------------
cci <- slopes[sapply(slopes, function(x) x$info$mask == "CCI")]

for (var in c("LAI","FPAR")){
  vv <- cci[sapply(cci, function(x) x$info$var == var)]
  if (length(vv) < 2) next

  for (i in seq_along(vv)){
    for (j in seq_along(vv)){
      if (i >= j) next

      d  <- vv[[i]]$r - vv[[j]]$r
      of <- file.path(
        OUTP,
        sprintf("%s_CCI_%s_minus_%s.png",
                var, vv[[i]]$info$tau, vv[[j]]$info$tau)
      )

      ttl <- sprintf("%s CCI: τ=%s − τ=%s",
                     var, vv[[i]]$info$tau, vv[[j]]$info$tau)

      plot_diff(d, ttl, of)
    }
  }
}

# ---------------------------------------------------------------------
# GLC − CCI
# ---------------------------------------------------------------------
glc <- slopes[sapply(slopes, function(x) x$info$mask == "GLC")]

for (var in c("LAI","FPAR")){
  gvar   <- glc[sapply(glc, function(x) x$info$var == var)]
  ccivar <- cci[sapply(cci, function(x) x$info$var == var)]

  for (g in gvar){
    t <- g$info$tau
    m <- ccivar[sapply(ccivar, function(x) x$info$tau == t)]
    if (length(m) == 0) next

    d  <- g$r - m[[1]]$r
    of <- file.path(OUTP, sprintf("%s_GLC_minus_CCI_tau_%s.png", var, t))

    ttl <- sprintf("%s: GLC − CCI (τ=%s)", var, t)

    plot_diff(d, ttl, of)
  }
}

cat("Difference plots saved in:", OUTP, "\n")
