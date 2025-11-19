## =============================================================================
# viz.R — Visualization and quicklook utilities for LAI/FPAR and mask products
#
# Purpose
#   Standardized plotting utilities for QA and reporting.
#
# Dependencies: terra, maps (optional), here, utils.R, geom.R
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R/utils.R"))
source(here("R/geom.R"))

# ==============================================================================
# Palettes
# ==============================================================================

pal_green <- function(n = 64) hcl.colors(n, "Greens", rev = TRUE)
pal_mask  <- c("#f0f0f0", "#d73027")   # keep / drop
col_na    <- "#bdbdbd"


# ==============================================================================
# Directory helper
# ==============================================================================

viz_dir_for <- function(cfg, step_tag) {
  d <- here(cfg$paths$viz_dir, step_tag)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}


# ==============================================================================
# Grid and coast overlays
# ==============================================================================

.add_overlays <- function(r,
                          grid_step = 30,
                          grid_col = "grey85",
                          coast_col = "grey25",
                          lwd = 0.6) {

  abline(h = seq(-90, 90, grid_step),
         v = seq(-180, 180, grid_step),
         col = grid_col, lwd = lwd)

  if (.have_maps) {
    ex <- ext(r)
    mp <- maps::map(
      "world",
      xlim = c(ex[1], ex[2]),
      ylim = c(ex[3], ex[4]),
      plot = FALSE
    )
    lines(mp$x, mp$y, col = coast_col, lwd = lwd)
  }
}


# ==============================================================================
# Generic raster plotters
# ==============================================================================

viz_png <- function(r,
                    outfile,
                    title = NULL,
                    zlim = NULL,
                    maxcell = 2e6,
                    categorical = FALSE) {

  if (file.exists(outfile)) return(invisible(outfile))

  png(outfile, width = 1600, height = 900, res = 150)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)

  par(mar = c(3, 3, 3, 5))
  if (categorical) r <- as.factor(r)

  terra::plot(
    r,
    main = title %||% basename(outfile),
    zlim = zlim,
    maxcell = maxcell,
    axes = FALSE,
    box = FALSE
  )
  box()
  invisible(outfile)
}


viz_hist <- function(r,
                     outfile,
                     breaks = NULL,
                     maxsamp = 1e6,
                     main = NULL,
                     xlab = "value") {

  if (file.exists(outfile)) return(invisible(outfile))

  png(outfile, width = 1400, height = 900, res = 150)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  par(mar = c(4, 4, 3, 1))

  rr <- try(terra::global(r, "range", na.rm = TRUE), silent = TRUE)
  if (inherits(rr, "try-error") || any(!is.finite(rr))) {
    plot.new()
    title(main = paste0(main %||% basename(outfile), " (no data)"))
    return(invisible(outfile))
  }

  vals <- unlist(rr[1, ], use.names = FALSE)
  lo   <- vals[1]; hi <- vals[2]
  if (hi <= lo) hi <- lo + max(1e-8, 1e-6 * abs(lo))

  breaks <- if (is.null(breaks)) pretty(c(lo, hi), 40) else {
    eps <- max(1e-8, 1e-6 * abs(hi - lo))
    c(min(breaks, lo - eps), max(breaks, hi + eps))
  }

  hh <- terra::hist(r, plot = FALSE, breaks = breaks, maxsamp = maxsamp, na.rm = TRUE)
  if (is.null(hh$counts) || !length(hh$counts)) {
    plot.new(); title(main = paste0(main %||% basename(outfile), " (empty)"))
    return(invisible(outfile))
  }

  plot(hh,
       main = main %||% basename(outfile),
       xlab = xlab)
  invisible(outfile)
}


# ==============================================================================
# AOI utilities
# ==============================================================================

aoi_extents <- function(cfg, drop_global = FALSE) {
  stopifnot(!is.null(cfg$aois))
  exts <- lapply(cfg$aois, function(a)
    ext(a$lon_min, a$lon_max, a$lat_min, a$lat_max)
  )
  if (drop_global) exts <- exts[setdiff(names(exts), "global")]
  exts
}


# ==============================================================================
# Mask visualization (global or AOI)
# ==============================================================================

plot_mask <- function(r, title, file, pal = pal_mask) {
  r2 <- terra::ifel(r >= 1, 1, 0)

  png(file, 1400, 700, res = 120)
  op <- par(mar = c(3, 3, 3, 6))
  on.exit({ par(op); dev.off() }, add = TRUE)

  terra::plot(r2, main = title,
              col = pal, breaks = c(-0.5, 0.5, 1.5),
              legend = FALSE, axes = TRUE, box = TRUE)

  legend("bottomleft", fill = pal, legend = c("0 keep", "1 drop"), bty = "n")
}


# ==============================================================================
# Before/after LAI/FPAR quicklooks
# ==============================================================================

quicklook_before_after <- function(rb, ra, ym, title, ql_dir, zlim, down = 4L) {

  dir.create(ql_dir, TRUE, showWarnings = FALSE)

  rb <- if (down > 1L) terra::aggregate(rb, down, mean, na.rm = TRUE) else rb
  ra <- if (down > 1L) terra::aggregate(ra, down, mean, na.rm = TRUE) else ra

  png(file.path(ql_dir, sprintf("quicklook_%s_masked_%s.png", title, ym)),
      width = 1400, height = 700, res = 120)

  op <- par(mfrow = c(1, 2), oma = c(0, 0, 2.2, 0), mar = c(3, 3, 3, 6))
  on.exit({ par(op); dev.off() }, add = TRUE)

  terra::plot(rb, main = sprintf("%s %s (before)", title, ym),
              col = pal_green(64), colNA = col_na, zlim = zlim,
              axes = TRUE, legend = TRUE)
  .add_overlays(rb)

  terra::plot(ra, main = sprintf("%s %s (masked)", title, ym),
              col = pal_green(64), colNA = col_na, zlim = zlim,
              axes = TRUE, legend = TRUE)
  .add_overlays(ra)

  mtext("Masked 0.05° quicklook", side = 3, outer = TRUE, cex = 1.05)
}


quicklook_after_full <- function(ra, ym, title, ql_dir, zlim, down = 4L) {

  dir.create(ql_dir, TRUE, showWarnings = FALSE)
  rr <- if (down > 1L) terra::aggregate(ra, down, mean, na.rm = TRUE) else ra

  png(file.path(ql_dir, sprintf("quicklook_%s_masked_full_%s.png", title, ym)),
      width = 1400, height = 700, res = 120)

  op <- par(oma = c(0, 0, 2, 0), mar = c(3, 3, 3, 6))
  on.exit({ par(op); dev.off() }, add = TRUE)

  terra::plot(rr, main = sprintf("%s %s (masked)", title, ym),
              col = pal_green(64), colNA = col_na, zlim = zlim,
              axes = TRUE, legend = TRUE)
  .add_overlays(rr)
}


# ==============================================================================
# Fractional cover quicklooks (cropland/urban)
# ==============================================================================

ql_write_two_panels <- function(r,
                                year,
                                title,
                                out_png,
                                zlim = c(0, 1),
                                pal = pal_green(64),
                                colNA = col_na) {

  .have_maps <- requireNamespace("maps", quietly = TRUE)

  add_overlay <- function() {
    abline(h = seq(-90, 90, 30),
           v = seq(-180, 180, 30),
           col = "grey85", lwd = 0.6)
    if (.have_maps) {
      mp <- maps::map("world", plot = FALSE)
      lines(mp$x, mp$y, col = "grey25", lwd = 0.6)
    }
  }

  png(out_png, width = 1400, height = 700, res = 120)
  op <- par(mfrow = c(1, 2), oma = c(2.2, 2.2, 3, 5), mar = c(3, 3, 2.5, 6))
  on.exit({ par(op); dev.off() }, add = TRUE)

  terra::plot(r[["frac_cropland"]], col = pal, colNA = colNA, zlim = zlim,
              axes = TRUE, legend = TRUE,
              plg = list(title = "Fraction [0,1]", cex = 0.8),
              main = sprintf("Cropland — %s %d", title, year))
  add_overlay()

  terra::plot(r[["frac_urban"]], col = pal, colNA = colNA, zlim = zlim,
              axes = TRUE, legend = TRUE,
              plg = list(title = "Fraction [0,1]", cex = 0.8),
              main = sprintf("Urban — %s %d", title, year))
  add_overlay()

  mtext("ESA-CCI/C3S 0.05° fractional cover quicklook",
        side = 3, outer = TRUE, cex = 1.05)
}


quicklook_all_aois <- function(frac,
                               year,
                               cfg,
                               ql_root,
                               down = 4L,
                               include_global = TRUE,
                               drop_global_key = FALSE) {

  stopifnot(inherits(frac, "SpatRaster"))
  dir.create(ql_root, TRUE, showWarnings = FALSE)

  x <- if (down > 1L) terra::aggregate(frac, down, mean, na.rm = TRUE) else frac

  if (include_global) {
    d <- file.path(ql_root, "global")
    dir.create(d, TRUE, showWarnings = FALSE)
    ql_write_two_panels(x, year, "Global",
                        file.path(d, sprintf("quicklook_global_%d.png", year)))
  }

  exts <- aoi_extents(cfg, drop_global_key)
  for (nm in names(exts)) {
    d <- file.path(ql_root, nm); dir.create(d, TRUE)
    rr <- try(terra::crop(x, exts[[nm]]), silent = TRUE)
    if (!inherits(rr, "try-error")) {
      ql_write_two_panels(rr, year, nm,
                          file.path(d, sprintf("quicklook_%s_%d.png", nm, year)))
    }
  }
}


# ==============================================================================
# Mask quicklooks
# ==============================================================================

ql_write_mask_two <- function(r_global,
                              r_local,
                              title_global,
                              title_local,
                              out_png,
                              col = pal_mask) {

  png(out_png, width = 1200, height = 600, res = 120)
  op <- par(mfrow = c(1, 2), oma = c(0, 0, 2, 0), mar = c(2, 2, 2, 5))
  on.exit({ par(op); dev.off() }, add = TRUE)

  brks <- c(-0.5, 0.5, 1.5)

  terra::plot(r_global, col = col, breaks = brks,
              legend = FALSE, main = title_global,
              axes = FALSE, box = TRUE)
  .add_overlays(r_global)
  legend("bottomleft", fill = col, legend = c("0 keep", "1 drop"), bty = "n")

  terra::plot(r_local, col = col, breaks = brks,
              legend = FALSE, main = title_local,
              axes = FALSE, box = TRUE)
  .add_overlays(r_local)
  legend("bottomleft", fill = col, legend = c("0 keep", "1 drop"), bty = "n")
}


quicklook_mask_all_aois <- function(mask,
                                    title,
                                    tag,
                                    cfg,
                                    ql_root,
                                    down = 4L,
                                    include_global = TRUE,
                                    drop_global_key = FALSE) {

  dir.create(ql_root, TRUE, showWarnings = FALSE)

  # Downscale by majority vote if needed
  x <- if (down > 1L) {
    maj <- function(v, ...) as.integer(mean(as.integer(v), na.rm = TRUE) >= 0.5)
    terra::aggregate(terra::ifel(mask, 1, 0), down, maj)
  } else {
    terra::ifel(mask, 1, 0)
  }

  if (include_global) {
    gdir <- file.path(ql_root, "global")
    dir.create(gdir, TRUE)
    ql_write_mask_two(
      x, x,
      paste(title, "— Global"),
      paste(title, "— Global"),
      file.path(gdir, sprintf("quicklook_mask_global_%s.png", tag))
    )
  }

  exts <- aoi_extents(cfg, drop_global_key)
  for (nm in names(exts)) {
    d <- file.path(ql_root, nm); dir.create(d, TRUE)
    rr <- try(terra::crop(x, exts[[nm]]), silent = TRUE)
    if (!inherits(rr, "try-error")) {
      ql_write_mask_two(
        x, rr,
        paste(title, "— Global"),
        paste(title, "—", nm),
        file.path(d, sprintf("quicklook_mask_%s_%s.png", nm, tag))
      )
    }
  }
}
