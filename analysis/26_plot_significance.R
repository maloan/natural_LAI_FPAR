#!/usr/bin/env Rscript
## =============================================================================
# 26_plot_significance.R — Mann–Kendall p-value maps
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(Kendall)
  library(here)
})

# ---------------------------------------------------------------------
# Root + output directory
# ---------------------------------------------------------------------
ROOT <- here()
OUTP <- here("analysis", "trend_plots", "significance")
dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)

pal <- hcl.colors(100, "batlow")

# ---------------------------------------------------------------------
# find monthly stacks (0.25°)
# ---------------------------------------------------------------------
files <- list.files(
  ROOT,
  pattern = "monthly_0p25_time\\.nc$",
  recursive = TRUE,
  full.names = TRUE
)

# ---------------------------------------------------------------------
# Mann–Kendall p-value raster
# ---------------------------------------------------------------------
mk_p <- function(ncfile){
  x <- rast(ncfile)        # layers = time
  vals <- values(x)

  p <- apply(vals, 1, function(ts){
    if (all(is.na(ts))) return(NA)
    r <- try(MannKendall(ts), silent = TRUE)
    if (inherits(r, "try-error")) return(NA)
    r$sl
  })

  r <- x[[1]]
  values(r) <- p
  r
}

# ---------------------------------------------------------------------
# main loop
# ---------------------------------------------------------------------
for (f in files){
  pmap <- mk_p(f)

  var  <- if (grepl("FPAR", f)) "FPAR" else "LAI"
  mask <- if (grepl("CCI", f)) "CCI" else if (grepl("GLC", f)) "GLC" else "UNMASKED"

  tau <- regmatches(f, regexpr("tau_0\\.[0-9]+", f))
  if (length(tau) == 0) tau <- "none"

  outdir <- file.path(OUTP, var, mask, tau)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  of <- file.path(
    outdir,
    sprintf("%s_%s_%s_MK_pvalues.png", var, mask, tau)
  )

  png(of, width = 1800, height = 900, res = 150)
  plot(
    pmap,
    col  = pal,
    main = sprintf("%s %s %s Mann–Kendall p-values", var, mask, tau)
  )
  dev.off()
}

cat("Significance maps saved to:", OUTP, "\n")
