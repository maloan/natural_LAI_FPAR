## =============================================================================
# viz.R — Visualization and quicklook utilities for LAI/FPAR and mask products
#
# Purpose
#   Provide standardized plotting utilities for QA and reporting, including
#   color palettes, histogram summaries, and global/AOI quicklook generation.
#
# Functions
#   viz_dir_for(cfg, step_tag)             — Create/resolve visualization subdir
#   viz_png(r, outfile, ...)               — Generic raster plotting
#   viz_hist(r, outfile, ...)              — Histogram helper with safe ranges
#   aoi_extents(cfg, drop_global)          — Build AOI extent list from cfg$aois
#   quicklook_before_after(rb, ra, ...)    — 2-panel before/after (masked)
#   quicklook_after_full(ra, ...)          — Single-panel masked preview
#   ql_write_two_panels(r, year, ...)      — Cropland/Urban 2-panel quicklook
#   quicklook_all_aois(frac, year, ...)    — Global + AOI quicklooks (fractions)
#   ql_write_mask_two(r_global, r_local, ...) — Mask 2-panel (global vs AOI)
#   quicklook_mask_all_aois(mask, ...)     — Global + AOI quicklooks (masks)
#
# Inputs
#   - SpatRaster objects (fractional or binary)
#   - cfg$aois : AOI metadata with lon_min/lon_max/lat_min/lat_max
#   - ql_root  : Output directory for PNG quicklooks
#
# Outputs
#   - PNG figures under <ql_root>/<AOI>/quicklook_*.png
#
# Dependencies
#   Packages: terra; maps (optional for coastlines)
#
# Notes
#   - No reliance on global variables like ROOT/out_dir/ql_title.
#   - All if/else use explicit curly braces for clarity and robustness.
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
})
source(file.path(ROOT, "R", "utils.R"))
source(file.path(ROOT, "R", "geom.R"))

# ----- Palettes and small constants -------------------------------------------

pal_green <- function(n = 64) {
  hcl.colors(n, "Greens", rev = TRUE)
}
pal_mask   <- c("#f0f0f0", "#d73027")  # 0 keep, 1 drop
col_na     <- "#bdbdbd"
.have_maps <- requireNamespace("maps", quietly = TRUE)

# ----- Directory helpers -------------------------------------------------------

viz_dir_for <- function(cfg, step_tag) {
  d <- file.path(path.expand(cfg$paths$viz_dir), step_tag)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  return(d)
}

# ----- Map overlays (graticule + coastlines when available) -------------------

.add_overlays <- function(r,
                          grid_step = 30,
                          grid_col = "grey85",
                          coast_col = "grey25",
                          lwd = 0.6) {
  abline(
    h = seq(-90, 90, by = grid_step),
    v = seq(-180, 180, by = grid_step),
    col = grid_col,
    lwd = lwd
  )
  if (.have_maps) {
    ex <- terra::ext(r)
    mp <- maps::map(
      "world",
      xlim = c(ex[1], ex[2]),
      ylim = c(ex[3], ex[4]),
      plot = FALSE
    )
    lines(mp$x, mp$y, col = coast_col, lwd = lwd)
  }
}

# ----- Generic plotters --------------------------------------------------------

viz_png <- function(r,
                    outfile,
                    title = NULL,
                    zlim = NULL,
                    maxcell = 2e6,
                    categorical = FALSE) {
  if (file.exists(outfile)) {
    return(invisible(outfile))
  }
  png(outfile,
      width = 1600,
      height = 900,
      res = 150)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  par(mar = c(3, 3, 3, 5))
  if (categorical) {
    r <- as.factor(r)
  }
  terra::plot(
    r,
    main = if (!is.null(title)) {
      title
    } else {
      basename(outfile)
    },
    zlim = zlim,
    maxcell = maxcell,
    axes = FALSE,
    box  = FALSE
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
  if (file.exists(outfile)) {
    return(invisible(outfile))
  }
  png(outfile,
      width = 1400,
      height = 900,
      res = 150)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  par(mar = c(4, 4, 3, 1))

  rr <- try(terra::global(r, "range", na.rm = TRUE), silent = TRUE)
  if (inherits(rr, "try-error")) {
    plot.new()
    title(main = paste0(main %||% basename(outfile), " (no data)"))
    return(invisible(outfile))
  }
  vals <- unlist(rr[1, ], use.names = FALSE)
  if (length(vals) < 2 || any(!is.finite(vals))) {
    plot.new()
    title(main = paste0(main %||% basename(outfile), " (no finite data)"))
    return(invisible(outfile))
  }
  lo <- as.numeric(vals[1])
  hi <- as.numeric(vals[2])
  if (!is.finite(lo) || !is.finite(hi)) {
    plot.new()
    title(main = paste0(main %||% basename(outfile), " (no finite data)"))
    return(invisible(outfile))
  }
  if (hi <= lo) {
    eps <- max(1e-8, 1e-6 * abs(lo))
    hi <- lo + eps
  }
  if (!is.null(breaks)) {
    eps <- max(1e-8, 1e-6 * abs(hi - lo))
    if (min(breaks) > lo) {
      breaks <- c(lo - eps, breaks)
    }
    if (max(breaks) < hi) {
      breaks <- c(breaks, hi + eps)
    }
  } else {
    breaks <- pretty(c(lo, hi), n = 40)
  }

  hh <- terra::hist(
    r,
    plot = FALSE,
    breaks = breaks,
    maxsamp = maxsamp,
    na.rm = TRUE
  )
  if (is.null(hh$counts) || length(hh$counts) == 0) {
    plot.new()
    title(main = paste0(main %||% basename(outfile), " (empty sample)"))
    return(invisible(outfile))
  }
  plot(hh, main = if (!is.null(main)) {
    main
  } else {
    basename(outfile)
  }, xlab = xlab)
  invisible(outfile)
}

# ----- AOI utilities -----------------------------------------------------------

aoi_extents <- function(cfg, drop_global = FALSE) {
  stopifnot(!is.null(cfg$aois))
  exts <- lapply(cfg$aois, function(a) {
    terra::ext(a$lon_min, a$lon_max, a$lat_min, a$lat_max)
  })
  if (drop_global && "global" %in% names(exts)) {
    exts <- exts[setdiff(names(exts), "global")]
  }
  return(exts)
}

# ----- LAI/FPAR masked quicklooks ---------------------------------------------

quicklook_before_after <- function(rb, ra, ym, title, ql_dir, zlim, down = 4L) {
  dir.create(ql_dir, TRUE, showWarnings = FALSE)
  if (down > 1L) {
    rb <- terra::aggregate(rb, down, mean, na.rm = TRUE)
    ra <- terra::aggregate(ra, down, mean, na.rm = TRUE)
  }
  png(
    file.path(ql_dir, sprintf(
      "quicklook_%s_masked_%s.png", title, ym
    )),
    width = 1400,
    height = 700,
    res = 120
  )
  op <- par(
    mfrow = c(1, 2),
    oma = c(0, 0, 2.2, 0),
    mar = c(3, 3, 3, 6)
  )
  terra::plot(
    rb,
    main = sprintf("%s %s (before)", title, ym),
    col = pal_green(64),
    colNA = col_na,
    zlim = zlim,
    axes = TRUE,
    legend = TRUE
  )
  .add_overlays(rb)
  mtext("Longitude (°E)", side = 1, line = 2)
  mtext("Latitude (°N)", side = 2, line = 2)
  terra::plot(
    ra,
    main = sprintf("%s %s (masked)", title, ym),
    col = pal_green(64),
    colNA = col_na,
    zlim = zlim,
    axes = TRUE,
    legend = TRUE
  )
  .add_overlays(ra)
  mtext("Longitude (°E)", side = 1, line = 2)
  mtext("Latitude (°N)", side = 2, line = 2)
  mtext(
    "Masked 0.05° quicklook",
    side = 3,
    outer = TRUE,
    cex = 1.05
  )
  par(op)
  dev.off()
}

quicklook_after_full <- function(ra, ym, title, ql_dir, zlim, down = 4L) {
  dir.create(ql_dir, TRUE, showWarnings = FALSE)
  rr <- if (down > 1L) {
    terra::aggregate(ra, down, mean, na.rm = TRUE)
  } else {
    ra
  }
  png(
    file.path(
      ql_dir,
      sprintf("quicklook_%s_masked_full_%s.png", title, ym)
    ),
    width = 1400,
    height = 700,
    res = 120
  )
  op <- par(oma = c(0, 0, 2, 0), mar = c(3, 3, 3, 6))
  terra::plot(
    rr,
    main = sprintf("%s %s (masked)", title, ym),
    col = pal_green(64),
    colNA = col_na,
    zlim = zlim,
    axes = TRUE,
    legend = TRUE
  )
  .add_overlays(rr)
  mtext("Longitude (°E)", side = 1, line = 2)
  mtext("Latitude (°N)", side = 2, line = 2)
  par(op)
  dev.off()
}

# ----- Fractional cover quicklooks (cropland/urban) ----------------------------

ql_write_two_panels <- function(r,
                                year,
                                title,
                                out_png,
                                zlim = c(0, 1),
                                pal  = pal_green(64),
                                colNA = col_na) {
  .have_maps_local <- requireNamespace("maps", quietly = TRUE)
  add_overlays_local <- function() {
    abline(
      h = seq(-90, 90, by = 30),
      v = seq(-180, 180, by = 30),
      col = "grey85",
      lwd = 0.6
    )
    if (.have_maps_local) {
      mp <- maps::map("world", plot = FALSE)
      lines(mp$x, mp$y, col = "grey25", lwd = 0.6)
    }
  }
  png(out_png,
      width = 1400,
      height = 700,
      res = 120)
  op <- par(
    mfrow = c(1, 2),
    oma = c(2.2, 2.2, 3, 5),
    mar = c(3, 3, 2.5, 6)
  )
  terra::plot(
    r[["frac_cropland"]],
    col = pal,
    colNA = colNA,
    zlim = zlim,
    axes = TRUE,
    legend = TRUE,
    plg = list(title = "Fraction [0,1]", cex = 0.8),
    main = sprintf("Cropland — %s %d", title, year)
  )
  add_overlays_local()
  mtext("Longitude (°E)", side = 1, line = -1)
  mtext("Latitude (°N)", side = 2, line = 2)
  terra::plot(
    r[["frac_urban"]],
    col = pal,
    colNA = colNA,
    zlim = zlim,
    axes = TRUE,
    legend = TRUE,
    plg = list(title = "Fraction [0,1]", cex = 0.8),
    main = sprintf("Urban — %s %d", title, year)
  )
  add_overlays_local()
  mtext("Longitude (°E)", side = 1, line = -1)
  mtext("Latitude (°N)", side = 2, line = 2)
  mtext(
    "ESA-CCI/C3S → 0.05° fractional cover quicklook",
    side = 3,
    outer = TRUE,
    cex = 1.05
  )
  par(op)
  dev.off()
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
  x <- if (down > 1L) {
    terra::aggregate(frac,
                     fact = down,
                     fun = mean,
                     na.rm = TRUE)
  } else {
    frac
  }
  if (isTRUE(include_global)) {
    d <- file.path(ql_root, "global")
    dir.create(d, TRUE, showWarnings = FALSE)
    ql_write_two_panels(
      r = x,
      year = year,
      title = "Global",
      out_png = file.path(d, sprintf("quicklook_global_%d.png", year))
    )
  }
  exts <- aoi_extents(cfg, drop_global = drop_global_key)
  for (nm in names(exts)) {
    d <- file.path(ql_root, nm)
    dir.create(d, TRUE, showWarnings = FALSE)
    rr <- try(terra::crop(x, exts[[nm]]), silent = TRUE)
    if (inherits(rr, "try-error")) {
      next
    } else {
      ql_write_two_panels(
        r = rr,
        year = year,
        title = nm,
        out_png = file.path(d, sprintf("quicklook_%s_%d.png", nm, year))
      )
    }
  }
}

# ----- Mask quicklooks (0 keep = grey, 1 drop = red) ---------------------------

ql_write_mask_two <- function(r_global,
                              r_local,
                              title_global,
                              title_local,
                              out_png,
                              col = pal_mask) {
  brks <- c(-0.5, 0.5, 1.5)
  png(out_png,
      width = 1200,
      height = 600,
      res = 120)
  op <- par(
    mfrow = c(1, 2),
    oma = c(0, 0, 2, 0),
    mar = c(2, 2, 2, 5)
  )
  terra::plot(
    r_global,
    col = col,
    breaks = brks,
    legend = FALSE,
    main = title_global,
    axes = FALSE,
    box = TRUE
  )
  .add_overlays(r_global)
  legend(
    "bottomleft",
    fill = col,
    legend = c("0 keep", "1 drop"),
    bty = "n"
  )
  terra::plot(
    r_local,
    col = col,
    breaks = brks,
    legend = FALSE,
    main = title_local,
    axes = FALSE,
    box = TRUE
  )
  .add_overlays(r_local)
  legend(
    "bottomleft",
    fill = col,
    legend = c("0 keep", "1 drop"),
    bty = "n"
  )
  par(op)
  dev.off()
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
  x <- if (down > 1L) {
    maj <- function(v, ...) {
      as.integer(mean(as.integer(v), na.rm = TRUE) >= 0.5)
    }
    terra::aggregate(terra::ifel(mask, 1L, 0L), down, maj)
  } else {
    terra::ifel(mask, 1L, 0L)
  }

  if (isTRUE(include_global)) {
    d <- file.path(ql_root, "global")
    dir.create(d, TRUE, showWarnings = FALSE)
    ql_write_mask_two(
      r_global     = x,
      r_local      = x,
      title_global = paste(title, "— Global"),
      title_local  = paste(title, "— Global"),
      out_png      = file.path(d, sprintf("quicklook_mask_global_%s.png", tag))
    )
  }

  exts <- aoi_extents(cfg, drop_global = drop_global_key)
  for (nm in names(exts)) {
    d <- file.path(ql_root, nm)
    dir.create(d, TRUE, showWarnings = FALSE)
    rr <- try(terra::crop(x, exts[[nm]]), silent = TRUE)
    if (inherits(rr, "try-error")) {
      next
    } else {
      ql_write_mask_two(
        r_global     = x,
        r_local      = rr,
        title_global = paste(title, "— Global"),
        title_local  = paste(title, "—", nm),
        out_png      = file.path(d, sprintf(
          "quicklook_mask_%s_%s.png", nm, tag
        ))
      )
    }
  }
}
