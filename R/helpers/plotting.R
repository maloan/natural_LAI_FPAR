## =============================================================================
## plotting.R — Visualization and quicklook utilities for LAI/FPAR and mask products
## =============================================================================

# ------------------------------------------------------------------------------
# Palettes and shared style
# ------------------------------------------------------------------------------
pal_green <- function(n = 64) {
  hcl.colors(n, "Greens", rev = TRUE)
}
pal_mask <- c("#f0f0f0", "#d73027") # 0 keep, 1 drop
col_na <- "#bdbdbd"

# ------------------------------------------------------------------------------
# Internal overlays for base plots
# ------------------------------------------------------------------------------
.have_maps <- requireNamespace("maps", quietly = TRUE)

.add_graticule <- function() {
  abline(
    h = seq(-90, 90, by = 30),
    v = seq(-180, 180, by = 30),
    col = "grey85",
    lwd = 0.6
  )
}

.add_coastlines <- function() {
  if (.have_maps) {
    mp <- maps::map("world", plot = FALSE)
    lines(mp$x, mp$y, col = "grey25", lwd = 0.6)
  }
}

.add_overlays <- function(r = NULL) {
  # `r` kept for backward compatibility; not used.
  .add_graticule()
  .add_coastlines()
}

.add_axis_labels <- function(line_x = 2, line_y = 2) {
  mtext("Longitude (°E)", 1, line = line_x)
  mtext("Latitude (°N)", 2, line = line_y)
}


# ============================================================================== 
# AOI utilities
# ==============================================================================

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

# ============================================================================== 
# Mask plotting + quicklooks (base graphics)
# ==============================================================================

quicklook_before_after <- function(rb, ra, ym, title, ql_dir, down = 2L) {
  stopifnot(inherits(rb, "SpatRaster"), inherits(ra, "SpatRaster"))
  dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

  if (down > 1L) {
    rb <- terra::aggregate(rb, down, mean, na.rm = TRUE)
    ra <- terra::aggregate(ra, down, mean, na.rm = TRUE)
  }

  out <- file.path(ql_dir, sprintf("quicklook_%s_masked_%s.png", title, ym))
  png(out, width = 1400, height = 700, res = 120)
  op <- par(mfrow = c(1, 2), oma = c(0, 0, 2.2, 0), mar = c(3, 3, 3, 6))
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
  )

  terra::plot(
    rb,
    main   = sprintf("%s %s (before)", title, ym),
    col    = pal_green(64),
    colNA  = col_na,
    axes   = TRUE,
    legend = TRUE
  )
  .add_overlays(rb)
  .add_axis_labels()

  terra::plot(
    ra,
    main   = sprintf("%s %s (masked)", title, ym),
    col    = pal_green(64),
    colNA  = col_na,
    axes   = TRUE,
    legend = TRUE
  )
  .add_overlays(ra)
  .add_axis_labels()

  mtext("Masked 0.05° quicklook", 3, outer = TRUE, cex = 1.05)
  invisible(out)
}

quicklook_after_full <- function(ra, ym, title, ql_dir, down = 2L) {
  stopifnot(inherits(ra, "SpatRaster"))
  dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

  rr <- if (down > 1L) terra::aggregate(ra, down, mean, na.rm = TRUE) else ra

  out <- file.path(
    ql_dir, sprintf("quicklook_%s_masked_full_%s.png", title, ym)
  )
  png(out, width = 1400, height = 700, res = 120)
  op <- par(oma = c(0, 0, 2, 0), mar = c(3, 3, 3, 6))
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
  )

  terra::plot(
    rr,
    main   = sprintf("%s %s (masked)", title, ym),
    col    = pal_green(64),
    colNA  = col_na,
    axes   = TRUE,
    legend = TRUE
  )
  .add_overlays(rr)
  .add_axis_labels()
  invisible(out)
}

# ============================================================================== 
# Fractional cover quicklooks (two-panel, base graphics)
# ==============================================================================

ql_write_two_panels <- function(r,
                                year,
                                title,
                                out_png,
                                zlim = c(0, 1),
                                pal = pal_green(64),
                                colNA = col_na) {
  stopifnot(inherits(r, "SpatRaster"))
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  stopifnot(all(c("frac_cropland", "frac_urban") %in% names(r)))

  png(out_png, width = 1400, height = 700, res = 120)
  op <- par(mfrow = c(1, 2), oma = c(2.2, 2.2, 3, 5), mar = c(3, 3, 2.5, 6))
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
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
  .add_overlays(r)
  .add_axis_labels(line_x = -1)

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
  .add_overlays(r)
  .add_axis_labels(line_x = -1)

  mtext(
    "ESA-CCI/C3S → 0.05° fractional cover quicklook", 3,
    outer = TRUE, cex = 1.05
  )
  invisible(out_png)
}

quicklook_all_aois <- function(frac,
                               year,
                               cfg,
                               ql_root,
                               down = 2L,
                               include_global = TRUE,
                               drop_global_key = FALSE) {
  stopifnot(inherits(frac, "SpatRaster"))
  dir.create(ql_root, recursive = TRUE, showWarnings = FALSE)

  x <- if (down > 1L) {
    terra::aggregate(frac, fact = down, fun = mean, na.rm = TRUE)
  } else {
    frac
  }

  if (isTRUE(include_global)) {
    d <- file.path(ql_root, "global")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    ql_write_two_panels(
      r = x, year = year, title = "Global",
      out_png = file.path(d, sprintf("quicklook_global_%d.png", year))
    )
  }

  exts <- aoi_extents(cfg, drop_global = drop_global_key)
  for (nm in names(exts)) {
    d <- file.path(ql_root, nm)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    rr <- try(terra::crop(x, exts[[nm]]), silent = TRUE)
    if (inherits(rr, "try-error")) next
    ql_write_two_panels(
      r = rr, year = year, title = nm,
      out_png = file.path(d, sprintf("quicklook_%s_%d.png", nm, year))
    )
  }

  invisible(ql_root)
}

# ============================================================================== 
# Mask quicklooks across AOIs
# ==============================================================================

ql_write_mask_two <- function(r_global,
                              r_local,
                              title_global,
                              title_local,
                              out_png,
                              col = pal_mask) {
  stopifnot(inherits(r_global, "SpatRaster"), inherits(r_local, "SpatRaster"))
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  brks <- c(-0.5, 0.5, 1.5)
  png(out_png, width = 1200, height = 600, res = 120)
  op <- par(mfrow = c(1, 2), oma = c(0, 0, 2, 0), mar = c(2, 2, 2, 5))
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
  )

  terra::plot(r_global,
    col = col, breaks = brks, legend = FALSE,
    main = title_global, axes = FALSE, box = TRUE
  )
  .add_overlays(r_global)
  legend("bottomleft", fill = col, legend = c("0 keep", "1 drop"), bty = "n")

  terra::plot(r_local,
    col = col, breaks = brks, legend = FALSE,
    main = title_local, axes = FALSE, box = TRUE
  )
  .add_overlays(r_local)
  legend("bottomleft", fill = col, legend = c("0 keep", "1 drop"), bty = "n")

  invisible(out_png)
}

quicklook_mask_all_aois <- function(mask,
                                    title,
                                    tag,
                                    cfg,
                                    ql_root,
                                    down = 2L,
                                    include_global = TRUE,
                                    drop_global_key = FALSE) {
  stopifnot(inherits(mask, "SpatRaster"))
  dir.create(ql_root, recursive = TRUE, showWarnings = FALSE)

  # Ensure 0/1/NA using terra-native conditionals (no base ifelse on rasters)
  m01 <- terra::ifel(
    terra::is.na(mask), NA,
    terra::ifel(mask >= 1, 1, 0)
  )
  m01 <- terra::as.int(m01)

  x <- if (down > 1L) {
    maj <- function(v, ...) {
      if (all(is.na(v))) {
        return(NA_integer_)
      }
      as.integer(mean(v, na.rm = TRUE) >= 0.5)
    }
    terra::aggregate(m01, fact = down, fun = maj)
  } else {
    m01
  }

  if (isTRUE(include_global)) {
    gdir <- file.path(ql_root, "global")
    dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
    ql_write_mask_two(
      x, x,
      paste(title, "— Global"),
      paste(title, "— Global"),
      file.path(gdir, sprintf("quicklook_mask_global_%s.png", tag))
    )
  }

  exts <- aoi_extents(cfg, drop_global = drop_global_key)
  for (nm in names(exts)) {
    d <- file.path(ql_root, nm)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
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

  invisible(ql_root)
}

# ------------------------------------------------------------------------------
# write_quicklook_raster() — standard quicklook PNG for a SpatRaster
# ------------------------------------------------------------------------------
write_quicklook_raster <- function(r,
                                   out_png,
                                   title = NULL,
                                   width = 1100,
                                   height = 550,
                                   res = 96,
                                   palette = "Greens",
                                   ncol = 64,
                                   rev = TRUE,
                                   colNA = col_na,
                                   axes = TRUE,
                                   ...) {
  stopifnot(inherits(r, "SpatRaster"))
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  png(out_png, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)

  terra::plot(
    r,
    main  = title,
    col   = hcl.colors(ncol, palette = palette, rev = rev),
    colNA = colNA,
    axes  = axes,
    ...
  )
  .add_overlays(r)
  invisible(out_png)
}

# ============================================================================== 
# ggplot helpers
# ==============================================================================

.plot_require_gg <- function() {
  req <- c("ggplot2", "scales", "scico", "sf", "rnaturalearth")
  ok <- vapply(req, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop("ggplot helpers need packages: ", paste(req[!ok], collapse = ", "))
  }
  invisible(TRUE)
}

theme_pub <- function(base_size = 12, include_legend = FALSE, include_strip = FALSE) {
  .plot_require_gg()
  theme_list <- list(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = base_size + 1, face = "bold"),
    plot.subtitle    = element_text(size = 10),
    axis.title       = element_text(size = base_size - 1),
    axis.text        = element_text(size = 9)
  )
  if (include_strip) {
    theme_list$strip.text <- element_text(size = 10, face = "bold")
  }
  if (include_legend) {
    theme_list$legend.position <- "bottom"
    theme_list$legend.box <- "vertical"
    theme_list$legend.text <- element_text(size = 9)
  }
  theme_bw(base_size = base_size) + theme(!!!theme_list)
}

lab_deg <- function() {
  .plot_require_gg()
  scales::label_number(suffix = "°")
}
