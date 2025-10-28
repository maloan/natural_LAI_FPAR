suppressPackageStartupMessages({
  library(terra)
})

source(file.path(ROOT, "R", "00_utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))

pal_grass   <- hcl.colors(64, "Greens",  rev = TRUE)
pal_pasture <- hcl.colors(64, "YlOrBr",  rev = TRUE)
viz_dir_for <- function(CFG, step_tag) {
  d <- file.path(path.expand(CFG$paths$viz_dir), step_tag)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.add_overlays <- function(r) {
  # light graticule
  abline(h = seq(-90,  90, by = 30), v = seq(-180, 180, by = 30),
         col = "grey85", lwd = 0.6)
  # country outlines (maps fallback)
  if (.have_maps) {
    ex <- terra::ext(r)
    mp <- maps::map("world",
                    xlim = c(ex[1], ex[2]), ylim = c(ex[3], ex[4]),
                    plot = FALSE)
    lines(mp$x, mp$y, col = "grey25", lwd = 0.6)
  }
}
viz_png <- function(r,
                    outfile,
                    title = NULL,
                    zlim = NULL,
                    maxcell = 2e6,
                    categorical = FALSE) {
  if (file.exists(outfile))
    return(invisible(outfile))
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
  if (categorical)
    r <- as.factor(r)
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
  if (file.exists(outfile))
    return(invisible(outfile))
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

  # safer: get min/max as numerics
  rr <- try(terra::global(r, "range", na.rm = TRUE), silent = TRUE)
  # cols: min, max
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
    # constant raster; make a tiny span so hist works
    eps <- max(1e-8, 1e-6 * abs(lo))
    hi <- lo + eps
  }

  # ensure breaks span the actual data range
  if (!is.null(breaks)) {
    eps <- max(1e-8, 1e-6 * abs(hi - lo))
    if (min(breaks) > lo)
      breaks <- c(lo - eps, breaks)
    if (max(breaks) < hi)
      breaks <- c(breaks, hi + eps)
  } else {
    breaks <- pretty(c(lo, hi), n = 40)
  }

  # compute with terra (sampling), plot with base to avoid forwarding 'maxsamp' to graphics
  hh <- terra::hist(
    r,
    plot = FALSE,
    breaks = breaks,
    maxsamp = maxsamp,
    na.rm = TRUE
  )
  # If hh is empty for some reason, fallback message
  if (is.null(hh$counts) || length(hh$counts) == 0) {
    plot.new()
    title(main = paste0(main %||% basename(outfile), " (empty sample)"))
    return(invisible(outfile))
  }
  plot(hh, main = main %||% basename(outfile), xlab = xlab)
  invisible(outfile)
}

# --- Quicklook helpers (global + all AOIs) ------------------------------------

# Small, consistent palette
pal_green <- function(n = 64) hcl.colors(n, "Greens", rev = TRUE)
col_na    <- "#bdbdbd"

# Return a *named* list of SpatExtent for every AOI in cfg$aois.
# Set drop_global=TRUE to skip the 'global' AOI key if present.
aoi_extents <- function(cfg, drop_global = FALSE) {
  stopifnot(!is.null(cfg$aois))
  exts <- lapply(cfg$aois, function(a) {
    terra::ext(a$lon_min, a$lon_max, a$lat_min, a$lat_max)
  })
  if (drop_global && "global" %in% names(exts)) {
    exts <- exts[setdiff(names(exts), "global")]
  }
  exts
}

# Draw one 2-panel quicklook (cropland + urban) to a single PNG
# 'r' must have layers named 'frac_cropland' and 'frac_urban'
ql_write_two_panels <- function(r,
                                year,
                                title,
                                out_png,
                                zlim = c(0, 1),
                                pal  = pal_green(64),
                                colNA = col_na) {

  # keep this base-graphics only + 'maps' fallback (no sf/rnaturalearth deps)
  .have_maps <- requireNamespace("maps", quietly = TRUE)

  add_overlays <- function() {
    # light graticule
    abline(h = seq(-90,  90, by = 30), v = seq(-180, 180, by = 30),
           col = "grey85", lwd = 0.6)
    # coastlines
    if (.have_maps) {
      mp <- maps::map("world", plot = FALSE)
      lines(mp$x, mp$y, col = "grey25", lwd = 0.6)
    }
  }

  png(out_png, width = 1400, height = 700, res = 120)
  op <- par(mfrow = c(1, 2),
            oma = c(2.2, 2.2, 3, 5),
            mar = c(3, 3, 2.5, 6))

  terra::plot(r[["frac_cropland"]],
              col   = pal, colNA = colNA, zlim = zlim,
              axes  = TRUE, legend = TRUE,
              plg   = list(title = "Fraction [0,1]", cex = 0.8),
              main  = sprintf("Cropland — %s %d", title, year))
  add_overlays()
  mtext("Longitude (°E)", side = 1, line = -1)
  mtext("Latitude (°N)",  side = 2, line = 2)

  terra::plot(r[["frac_urban"]],
              col   = pal, colNA = colNA, zlim = zlim,
              axes  = TRUE, legend = TRUE,
              plg   = list(title = "Fraction [0,1]", cex = 0.8),
              main  = sprintf("Urban — %s %d", title, year))
  add_overlays()
  mtext("Longitude (°E)", side = 1, line = -1)
  mtext("Latitude (°N)",  side = 2, line = 2)

  mtext("ESA-CCI/C3S → 0.05° fractional cover quicklook",
        side = 3, outer = TRUE, cex = 1.05)
  par(op); dev.off()
}

# Main entry: write Global + per-AOI quicklooks for a fraction stack
# frac: SpatRaster with at least 'frac_cropland' and 'frac_urban'
# ql_root: directory where 'global/' and '<AOI_NAME>/' folders will be created
quicklook_all_aois <- function(frac,
                               year,
                               cfg,
                               ql_root,
                               down = 4L,
                               include_global = TRUE,
                               drop_global_key = FALSE) {
  stopifnot(inherits(frac, "SpatRaster"))
  dir.create(ql_root, TRUE, showWarnings = FALSE)

  # Optional downsample for speed
  x <- if (down > 1L) terra::aggregate(frac, fact = down, fun = mean, na.rm = TRUE) else frac

  # 1) Global
  if (isTRUE(include_global)) {
    d <- file.path(ql_root, "global")
    dir.create(d, TRUE, showWarnings = FALSE)
    ql_write_two_panels(
      r = x, year = year, title = "Global",
      out_png = file.path(d, sprintf("quicklook_global_%d.png", year))
    )
  }

  # 2) All AOIs from cfg$aois
  exts <- aoi_extents(cfg, drop_global = drop_global_key)
  for (nm in names(exts)) {
    d <- file.path(ql_root, nm)
    dir.create(d, TRUE, showWarnings = FALSE)
    rr <- try(terra::crop(x, exts[[nm]]), silent = TRUE)
    if (inherits(rr, "try-error")) next
    ql_write_two_panels(
      r = rr, year = year, title = nm,
      out_png = file.path(d, sprintf("quicklook_%s_%d.png", nm, year))
    )
  }
}

# --- Mask quicklooks (0 keep = grey, 1 drop = red) ----------------------------

pal_mask   <- c("#f0f0f0", "#d73027")  # 0 keep, 1 drop
.have_maps <- requireNamespace("maps", quietly = TRUE)

.add_outline <- function(r) {
  # light graticule
  abline(h = seq(-90,  90, by = 30), v = seq(-180, 180, by = 30),
         col = "grey85", lwd = 0.6)
  # country outlines (fallback via 'maps')
  if (.have_maps) {
    ex <- terra::ext(r)
    mp <- maps::map("world",
                    xlim = c(ex[1], ex[2]),
                    ylim = c(ex[3], ex[4]),
                    plot = FALSE)
    lines(mp$x, mp$y, col = "grey25", lwd = 0.6)
  }
}

# Draw a two-panel PNG from *two* rasters with different extents
# r_global: full-extent mask (0/1)
# r_local : cropped-to-AOI mask (0/1)
ql_write_mask_two <- function(r_global, r_local, title_global, title_local, out_png,
                              col = pal_mask) {
  brks <- c(-0.5, 0.5, 1.5)
  png(out_png, width = 1200, height = 600, res = 120)
  op <- par(mfrow = c(1, 2), oma = c(0,0,2,0), mar = c(2,2,2,5))

  terra::plot(r_global, col = col, breaks = brks, legend = FALSE,
              main = title_global, axes = FALSE, box = TRUE)
  .add_outline(r_global)
  legend("bottomleft", fill = col, legend = c("0 keep", "1 drop"), bty = "n")

  terra::plot(r_local, col = col, breaks = brks, legend = FALSE,
              main = title_local, axes = FALSE, box = TRUE)
  .add_outline(r_local)
  legend("bottomleft", fill = col, legend = c("0 keep", "1 drop"), bty = "n")

  par(op); dev.off()
}

# Global + all AOIs driver
# mask: logical/byte SpatRaster (1=drop, 0=keep)
quicklook_mask_all_aois <- function(mask, title, tag, cfg, ql_root, down = 4L,
                                    include_global = TRUE, drop_global_key = FALSE) {
  dir.create(ql_root, TRUE, showWarnings = FALSE)

  # aggregate (majority) if requested
  x <- if (down > 1L) {
    maj <- function(v, ...) as.integer(mean(as.integer(v), na.rm = TRUE) >= 0.5)
    terra::aggregate(terra::ifel(mask, 1L, 0L), down, maj)
  } else {
    terra::ifel(mask, 1L, 0L)
  }

  # 1) Global
  if (isTRUE(include_global)) {
    d <- file.path(ql_root, "global"); dir.create(d, TRUE, showWarnings = FALSE)
    ql_write_mask_two(
      r_global     = x,
      r_local      = x,
      title_global = paste(title, "— Global"),
      title_local  = paste(title, "— Global"),
      out_png      = file.path(d, sprintf("quicklook_mask_global_%s.png", tag))
    )
  }

  # 2) All AOIs
  exts <- aoi_extents(cfg, drop_global = drop_global_key)
  for (nm in names(exts)) {
    d <- file.path(ql_root, nm); dir.create(d, TRUE, showWarnings = FALSE)
    rr <- try(terra::crop(x, exts[[nm]]), silent = TRUE)
    if (inherits(rr, "try-error")) next

    ql_write_mask_two(
      r_global     = x,
      r_local      = rr,
      title_global = paste(title, "— Global"),
      title_local  = paste(title, "—", nm),
      out_png      = file.path(d, sprintf("quicklook_mask_%s_%s.png", nm, tag))
    )
  }
}


