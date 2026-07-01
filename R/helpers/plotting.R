# =============================================================================
# plotting.R — Helper functions for plotting and visualizing raster data,
# trends, and quicklooks using base graphics and ggplot2.
# =============================================================================

# ------------------------------------------------------------------------------
# Palettes and shared style
# ------------------------------------------------------------------------------
pal_green <- function(n = 64) {
  # Return a green color palette for plotting, using HCL colors.
  hcl.colors(n, "Greens", rev = TRUE)
}
pal_mask <- c("#f0f0f0", "#d73027") # 0 keep, 1 drop
col_na <- "#bdbdbd"
lat_labels <- function(x) {
  # Return latitude labels in degrees North/South format for plotting.
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N")))
}


.plot_require_gg <- function() {
  # Check if required ggplot2-related packages are available, and stop with an
  # error message if any are missing.
  req <- c("ggplot2", "scales", "scico", "sf", "rnaturalearth")
  ok <- vapply(req, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop("ggplot helpers need packages: ", paste(req[!ok], collapse = ", "))
  }
  invisible(TRUE)
}

theme_pub <- function(base_size = 12,
                      include_legend = FALSE,
                      include_strip = FALSE) {
  # Return a ggplot2 theme, with options to include legend and strip text. The
  # base font size can be adjusted.
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
theme_map <- function(base_size = 10,
                      legend_position = NULL) {
  # Return a ggplot2 theme suitable for map plots, with options to set the base
  # font size and legend position. The theme includes grid lines, axis text,
  # plot title, and legend styling.
  theme_list <- list(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.5, "cm")
  )
  if (!is.null(legend_position)) {
    theme_list$legend.position <- legend_position
  }
  theme_bw(base_size = base_size) + theme(!!!theme_list)
}
lon_breaks_60 <- seq(-180, 180, by = 60)
lat_breaks_crop <- c(-60, -30, 0, 30, 60)

countries <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") |>
  sf::st_transform(4326)
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf") |>
  sf::st_transform(4326)
# ------------------------------------------------------------------------------
# Internal overlays for base plots
# ------------------------------------------------------------------------------
.have_maps <- requireNamespace("maps", quietly = TRUE)

geom_sig_marker <- function(data = NULL,
                            mapping = aes(x = x, y = y),
                            sig_col = "sig_flag",
                            y_offset = 0,
                            size = 3.5,
                            color = "black",
                            ...) {
  # Add significance markers to a ggplot2 plot, based on a specified
  # significance column (sig_col) in the data. The markers are placed at the top
  # of the plot (y = Inf) with an optional vertical offset (y_offset). The size
  # and color of the markers can be customized.
  geom_text(
    mapping = aes(label = .data[[!!sig_col]]),
    data = data,
    y = Inf,
    vjust = -0.5,
    size = size,
    color = color,
    inherit.aes = TRUE,
    ...
  )
}

lab_deg <- function() {
  # Return a label formatter for degrees, using the scales package.
  .plot_require_gg()
  scales::label_number(suffix = "°")
}
geom_sig_col <- function(aes_sig = aes(fill_alpha = sig),
                         alpha_sig = 1.0,
                         alpha_nonsig = 0.3,
                         ...) {
  # Add significance-based alpha transparency to a ggplot2 plot, using a
  # specified aesthetic mapping (aes_sig) for the significance column (sig). The
  # alpha values for significant and non-significant data points can be
  # customized. Additional parameters can be passed via '...'.
  scale_alpha_manual(
    values = c("TRUE" = alpha_sig, "FALSE" = alpha_nonsig),
    guide = "none",
    ...
  )
}


highlight_ci_crossing <- function(ci_lower,
                                  ci_upper,
                                  color_nonsig = "grey70",
                                  color_sig = NA,
                                  linewidth = 1.2) {
  # Return a color vector for plotting, highlighting areas where the confidence
  # interval (ci_lower, ci_upper) crosses zero. Non-significant areas are
  # colored with color_nonsig, while significant areas are colored with
  # color_sig. The linewidth parameter can be used for line plots.
  crosses_zero <- ci_lower * ci_upper <= 0
  ifelse(crosses_zero, color_nonsig, color_sig)
}


sym_q_lim <- function(r, q = 0.99) {
  # Calculate symmetric quantile limits for a raster (r) based on the specified
  # quantile (q). The function returns a vector of two values: the negative and
  # positive limits, which can be used for color scaling in plots.
  v <- values(r, mat = FALSE)
  v <- v[is.finite(v)]
  lim <- as.numeric(stats::quantile(abs(v), probs = q, na.rm = TRUE))
  c(-lim, lim)
}

to_df <- function(r, name = "z") {
  # Convert a SpatRaster (r) to a data frame with longitude, latitude, and a
  # specified value column name. The function retains NA values in the output
  # data frame.
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df) <- c("lon", "lat", name)
  df
}

lon_labels <- function(x) {
  #  Return longitude labels in degrees East/West format for plotting.
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°W"), paste0(x, "°E")))
}

row_label <- function(txt) {
  # Return a rotated text grob for use as a row label in patchwork plots.
  patchwork::wrap_elements(full = grid::textGrob(
    txt,
    rot = 90,
    gp = grid::gpar(fontface = "bold", fontsize = 11)
  ))
}

add_panel_tag <- function(p, tag) {
  # Add a panel tag (e.g., "(a)", "(b)") to a ggplot2 plot (p) using the
  # specified tag text. The tag is positioned in the top-left corner of the
  # plot.
  p +
    ggplot2::labs(tag = tag) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 11),
      plot.tag.position = c(0.02, 0.98)
    )
}


.add_graticule <- function() {
  # Add a graticule (latitude/longitude grid) to the current plot using base
  # graphics.
  abline(
    h = seq(-90, 90, by = 30),
    v = seq(-180, 180, by = 30),
    col = "grey85",
    lwd = 0.6
  )
}

.add_coastlines <- function() {
  # Add coastlines to the current plot using the 'maps' package, if available.
  if (.have_maps) {
    mp <- maps::map("world", plot = FALSE)
    lines(mp$x, mp$y, col = "grey25", lwd = 0.6)
  }
}

.add_overlays <- function() {
  # Add graticule and coastlines to the current plot.

  .add_graticule()
  .add_coastlines()
}

.add_axis_labels <- function(line_x = 2, line_y = 2) {
  # Add axis labels for longitude and latitude to the current plot using base
  # graphics.
  mtext("Longitude (°E)", 1, line = line_x)
  mtext("Latitude (°N)", 2, line = line_y)
}


# ------------------------------------------------------------------------------
# AOI utilities
# ------------------------------------------------------------------------------

aoi_extents <- function(cfg, drop_global = FALSE) {
  # Return a list of terra::ext objects for each AOI defined in the
  # configuration, optionally dropping the "global" AOI.
  stopifnot(!is.null(cfg$aois))
  exts <- lapply(cfg$aois, function(a) {
    terra::ext(a$lon_min, a$lon_max, a$lat_min, a$lat_max)
  })
  if (drop_global && "global" %in% names(exts)) {
    exts <- exts[setdiff(names(exts), "global")]
  }
  exts
}

# -------------------------------------------------------------------------------
# Mask plotting + quicklooks (base graphics)
# -------------------------------------------------------------------------------

quicklook_before_after <- function(rb, ra, ym, title, ql_dir, down = 2L) {
  # Generate a two-panel quicklook PNG comparing the original raster (rb) and
  # the masked raster (ra) for a given year-month (ym) and title, saving to
  # ql_dir. Optionally downsample by 'down' factor.
  stopifnot(inherits(rb, "SpatRaster"), inherits(ra, "SpatRaster"))
  dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

  if (down > 1L) {
    rb <- terra::aggregate(rb, down, mean, na.rm = TRUE)
    ra <- terra::aggregate(ra, down, mean, na.rm = TRUE)
  }

  out <- file.path(ql_dir, sprintf("quicklook_%s_masked_%s.png", title, ym))
  png(out,
    width = 1400,
    height = 700,
    res = 120
  )
  op <- par(
    mfrow = c(1, 2),
    oma = c(0, 0, 2.2, 0),
    mar = c(3, 3, 3, 6)
  )
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
  # Generate a quicklook PNG for the masked raster (ra) for a given year-month
  # (ym) and title, saving to ql_dir. Optionally downsample by 'down' factor.
  stopifnot(inherits(ra, "SpatRaster"))
  dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

  rr <- if (down > 1L) {
    terra::aggregate(ra, down, mean, na.rm = TRUE)
  } else {
    ra
  }

  out <- file.path(
    ql_dir,
    sprintf("quicklook_%s_masked_full_%s.png", title, ym)
  )
  png(out,
    width = 1400,
    height = 700,
    res = 120
  )
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

# ------------------------------------------------------------------------------
# Fractional cover quicklooks (two-panel, base graphics)
# ------------------------------------------------------------------------------

ql_write_two_panels <- function(r,
                                year,
                                title,
                                out_png,
                                zlim = c(0, 1),
                                pal = pal_green(64),
                                colNA = col_na) {
  # Generate a two-panel quicklook PNG for fractional cover rasters (cropland
  # and urban) for a given year and title, saving to out_png. Optionally specify
  # zlim, color palette, and NA color.
  stopifnot(inherits(r, "SpatRaster"))
  dir.create(dirname(out_png),
    recursive = TRUE,
    showWarnings = FALSE
  )
  stopifnot(all(c("frac_cropland", "frac_urban") %in% names(r)))

  png(out_png,
    width = 1400,
    height = 700,
    res = 120
  )
  op <- par(
    mfrow = c(1, 2),
    oma = c(2.2, 2.2, 3, 5),
    mar = c(3, 3, 2.5, 6)
  )
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
    "ESA-CCI/C3S → 0.05° fractional cover quicklook",
    3,
    outer = TRUE,
    cex = 1.05
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
  # Generate quicklook PNGs for fractional cover rasters across all AOIs defined
  # in the configuration, for a given year, saving to ql_root. Optionally
  # downsample by 'down' factor, include a global quicklook, and drop the
  # "global" AOI if specified.
  stopifnot(inherits(frac, "SpatRaster"))
  dir.create(ql_root, recursive = TRUE, showWarnings = FALSE)

  x <- if (down > 1L) {
    terra::aggregate(frac,
      fact = down,
      fun = mean,
      na.rm = TRUE
    )
  } else {
    frac
  }

  if (isTRUE(include_global)) {
    d <- file.path(ql_root, "global")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
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
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    rr <- try(terra::crop(x, exts[[nm]]), silent = TRUE)
    if (inherits(rr, "try-error")) {
      next
    }
    ql_write_two_panels(
      r = rr,
      year = year,
      title = nm,
      out_png = file.path(d, sprintf("quicklook_%s_%d.png", nm, year))
    )
  }

  invisible(ql_root)
}

# ------------------------------------------------------------------------------
# Mask quicklooks across AOIs
# ------------------------------------------------------------------------------

ql_write_mask_two <- function(r_global,
                              r_local,
                              title_global,
                              title_local,
                              out_png,
                              col = pal_mask) {
  # Generate a two-panel quicklook PNG comparing a global mask raster (r_global)
  # and a local AOI mask raster (r_local), with specified titles and output file
  # path. The color palette for the mask can be customized.
  stopifnot(
    inherits(r_global, "SpatRaster"),
    inherits(r_local, "SpatRaster")
  )
  dir.create(dirname(out_png),
    recursive = TRUE,
    showWarnings = FALSE
  )

  brks <- c(-0.5, 0.5, 1.5)
  png(out_png,
    width = 1200,
    height = 600,
    res = 120
  )
  op <- par(
    mfrow = c(1, 2),
    oma = c(0, 0, 2, 0),
    mar = c(2, 2, 2, 5)
  )
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
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
  # Generate quicklook PNGs for a mask raster across all AOIs defined in the
  # configuration, for a given title and tag, saving to ql_root. Optionally
  # downsample by 'down' factor, include a global quicklook, and drop the
  # "global" AOI if specified.
  stopifnot(inherits(mask, "SpatRaster"))
  dir.create(ql_root, recursive = TRUE, showWarnings = FALSE)

  # Ensure 0/1/NA using terra-native conditionals (no base ifelse on rasters)
  m01 <- terra::ifel(terra::is.na(mask), NA, terra::ifel(mask >= 1, 1, 0))
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
      x,
      x,
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
    if (inherits(rr, "try-error")) {
      next
    }

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
  # Write a quicklook PNG for a given raster (r) to the specified output path
  # (out_png), with optional title, dimensions, resolution, color palette,
  # number of colors, NA color, and axes display. Additional plotting parameters
  # can be passed via '...'.
  stopifnot(inherits(r, "SpatRaster"))
  dir.create(dirname(out_png),
    recursive = TRUE,
    showWarnings = FALSE
  )

  png(out_png,
    width = width,
    height = height,
    res = res
  )
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

# ------------------------------------------------------------------------------
# Plots for analysis scripts
# ------------------------------------------------------------------------------

plot_timeseries <- function(df, trend_df, annotation_df, theme_pub) {
  # Generate a time series plot using ggplot2, showing the original data (df),
  # trend lines (trend_df), and annotations (annotation_df). The plot is faceted
  # by metric and includes custom colors, line types, and labels. The theme_pub
  # function is used for consistent styling.
  cols <- c(
    "Unmasked" = "black",
    "CCI tau=0.05" = "#1b9e77",
    "CCI tau=0.1" = "#d95f02",
    "CCI tau=0.2" = "#7570b3",
    "GLC" = "#386cb0"
  )
  metric_labs <- c("Annual mean" = "(a) Annual mean", "Annual maximum" = "(b) Annual maximum")
  p <- ggplot() +
    geom_line(
      data = df,
      aes(year, value, colour = scenario, linetype = scenario),
      alpha = 0.4,
      linewidth = 0.4
    ) +
    geom_line(
      data = trend_df,
      aes(year, fit, colour = scenario, linetype = scenario),
      linewidth = 0.7,
      alpha = 1.0
    ) +
    geom_text(
      data = annotation_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 2.5,
      colour = "grey20",
      lineheight = 0.95
    ) +
    facet_wrap(
      ~metric,
      nrow = 1,
      strip.position = "top",
      scales = "free_y",
      labeller = labeller(metric = metric_labs)
    ) +
    scale_colour_manual(values = cols, name = "Mask Scenario") +
    scale_linetype_manual(
      values = c(
        "Unmasked" = "solid",
        "CCI tau=0.05" = "dashed",
        "CCI tau=0.1" = "dashed",
        "CCI tau=0.2" = "dashed",
        "GLC" = "dotdash"
      ),
      name = "Mask Scenario"
    ) +
    guides(colour = guide_legend(
      nrow = 1,
      override.aes = list(
        linetype = c("solid", "dashed", "dashed", "dashed", "dotdash"),
        alpha = 1
      )
    )) +
    labs(x = "Year", y = expression("Global mean LAI (" * m^2 ~ m^-2 * ")")) +
    theme_pub(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.spacing.x = unit(0.2, "lines"),
      legend.key.width = unit(1.5, "lines"),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5,
        color = "grey40"
      )
    )
  p
}


plot_seasonal_amplitude <- function(z_abs) {
  # Generate a zonal plot of mean absolute seasonal amplitude (z_abs) using
  # ggplot2, with latitude on the x-axis and mean amplitude on the y-axis. The
  # plot includes lines for different scenarios, custom colors, and labels. The
  # theme_pub function is used for consistent styling.
  col_abs <- c(
    "Unmasked" = "grey20",
    "CCI tau=0.05" = "#1b9e77",
    "CCI tau=0.1" = "#d95f02",
    "CCI tau=0.2" = "#7570b3",
    "GLC" = "#386cb0"
  )

  ylab_abs <- if (toupper(var) == "LAI") {
    expression("Mean absolute seasonal amplitude LAI (" * m^2 ~ m^-2 * ")")
  } else {
    "Mean absolute seasonal amplitude fAPAR (-)"
  }
  p <- ggplot(z_abs, aes(lat_band, mean_yearamp, colour = scenario)) +
    geom_hline(
      yintercept = 0,
      colour = "grey90",
      linewidth = 0.25
    ) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    scale_x_continuous(
      limits = c(-60, 90),
      breaks = seq(-90, 90, by = 30),
      labels = lat_labels
    ) +
    scale_colour_manual(values = col_abs) +
    labs(x = "Latitude", y = ylab_abs, colour = NULL) +
    theme_pub(include_legend = TRUE, include_strip = TRUE)
  p
}

plot_zonal_reltrend <- function(df_plot, ttl) {
  # Generate a zonal plot of relative trends (df_plot) using ggplot2, with
  # latitude on the x-axis and relative trend percentage per year on the y-axis.
  # The plot includes confidence bands, lines for different scenarios, custom
  # colors, and labels. The theme_pub function is used for consistent styling.
  lat_breaks_30 <- seq(-90, 90, by = 30)
  col_map <- c(
    "Unmasked" = "grey20",
    "CCI tau=0.05" = "#1b9e77",
    "CCI tau=0.1" = "#d95f02",
    "CCI tau=0.2" = "#7570b3",
    "GLC" = "#386cb0"
  )

  p <- ggplot(
    df_plot,
    aes(
      .data$lat_band,
      .data$reltrend_pct_per_year,
      colour = .data$scenario
    )
  ) +
    # Confidence band (shaded ribbon)
    geom_ribbon(
      aes(
        ymin = .data$ci_lower,
        ymax = .data$ci_upper,
        fill = .data$scenario
      ),
      alpha = 0.15,
      colour = NA,
      na.rm = TRUE
    ) +
    geom_hline(
      yintercept = 0,
      colour = "grey90",
      linewidth = 0.25
    ) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    scale_colour_manual(values = col_map, drop = FALSE) +
    scale_fill_manual(
      values = col_map,
      drop = FALSE,
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(-60, 90),
      breaks = lat_breaks_30,
      labels = lat_labels
    ) +
    labs(
      x = "Latitude",
      y = expression("Relative LAI trend (% yr"^
        {
          -1
        } * ")"),
      colour = NULL
    ) +
    theme_pub(include_legend = TRUE, include_strip = TRUE)
  p
}


plot_map <- function(df,
                     zcol,
                     lims,
                     title = NULL,
                     fill_title = NULL,
                     df_grey = NULL) {
  # Generate a map plot using ggplot2, displaying data from a data frame (df)
  # with longitude and latitude coordinates. The plot includes tiles colored by
  # the specified zcol variable, optional greyed-out areas (df_grey),
  # coastlines, and country borders. The color scale is set using the scico
  # package, and the plot is styled with a map theme.
  ggplot(df) +
    geom_tile(aes(
      x = .data$lon,
      y = .data$lat,
      fill = .data[[zcol]]
    )) +
    (if (!is.null(df_grey)) {
      geom_tile(
        data = df_grey,
        inherit.aes = FALSE,
        aes(x = .data$lon, y = .data$lat),
        fill = "grey70",
        alpha = 0.45
      )
    }) +
    geom_sf(
      data = coast,
      linewidth = 0.1,
      colour = "black",
      fill = NA,
      inherit.aes = FALSE
    ) +
    geom_sf(
      data = countries,
      fill = NA,
      color = "black",
      linewidth = 0.1
    ) +
    coord_sf(
      xlim = c(-180, 180),
      ylim = c(-60, 85),
      expand = FALSE
    ) +
    scale_x_continuous(breaks = seq(-180, 180, by = 60), labels = lon_labels) +
    scale_y_continuous(breaks = c(-60, -30, 0, 30, 60), labels = lat_labels) +
    scale_fill_scico(
      palette = "vik",
      direction = -1,
      limits = lims,
      midpoint = 0,
      oob = scales::squish,
      na.value = "grey98",
      name = fill_title
    ) +
    labs(title = title) +
    theme_map(legend_position = "bottom")
}

plot_zonal_diagnostics <- function(df, ycol, ttl, ylab, y_limits = NULL) {
  # Generate a zonal plot of diagnostics (df) using ggplot2, with latitude on
  # the x-axis and the specified ycol variable on the y-axis. The plot includes
  # lines for different scenarios, custom colors, confidence bands if available,
  # and labels. The theme_pub function is used for consistent styling.
  col_map <- c(
    "Unmasked" = "grey20",
    "CCI tau=0.05" = "#1b9e77",
    "CCI tau=0.1" = "#d95f02",
    "CCI tau=0.2" = "#7570b3",
    "GLC" = "#386cb0"
  )
  p <- ggplot(df, aes(.data$lat_band, .data[[ycol]], colour = .data$scenario)) +
    geom_hline(
      yintercept = 0,
      colour = "grey40",
      linewidth = 0.4
    ) +
    geom_line(linewidth = 0.5, na.rm = TRUE) +
    scale_colour_manual(values = col_map, drop = FALSE) +
    scale_x_continuous(
      limits = c(-60, 90),
      breaks = seq(-60, 90, by = 30),
      labels = lat_label_fn
    ) +
    labs(
      title = ttl,
      x = "Latitude",
      y = ylab,
      colour = "Mask Scenario"
    ) +
    theme_pub_fn(base_size = 10.5, include_legend = TRUE) +
    theme(
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5
      ),
      axis.text = element_text(size = 9.5),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "horizontal",
      legend.margin = margin(t = 10, b = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95", linewidth = 0.1)
    ) +
    guides(title = "Mask Scenario", colour = guide_legend(nrow = 5, byrow = TRUE))

  if ("ci_lower" %in% names(df)) {
    p <- p +
      geom_ribbon(
        aes(
          x = .data$lat_band,
          ymin = .data$ci_lower,
          ymax = .data$ci_upper,
          fill = .data$scenario,
          group = .data$scenario
        ),
        alpha = 0.08,
        linewidth = 0,
        na.rm = TRUE,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = col_map, drop = FALSE)
  }

  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }

  p
}

plot_lc_trend <- function(plot_tab, plot_long, scale_factor) {
  # Generate a land cover trend plot using ggplot2, showing trends and
  # confidence intervals for different land cover types. The plot includes
  # segments connecting unmasked and masked means, error bars for confidence
  # intervals, and points for trends. The x-axis is scaled based on the provided
  # scale_factor.
  finite_x <- c(plot_long$trend, plot_long$ci_lower, plot_long$ci_upper)
  finite_x <- finite_x[is.finite(finite_x)]
  x_min <- 0
  trend_xmax <- max(finite_x) + 0.3 * max(finite_x)
  area_label_x <- max(finite_x) + 0.3 * max(finite_x)
  area_header_x <- area_label_x
  x_range <- trend_xmax - x_min
  p <- ggplot() +
    geom_segment(
      data = plot_tab,
      aes(
        x = mean_unmasked,
        xend = mean_masked,
        y = lc_name,
        yend = lc_name
      ),
      linewidth = 0.7,
      colour = "grey45",
      lineend = "round"
    ) +
    geom_errorbar(
      data = plot_long,
      aes(xmin = ci_lower, xmax = ci_upper, y = lc_name),
      width = 0.12,
      colour = "grey50",
      alpha = 0.5,
      na.rm = TRUE
    ) +
    geom_point(
      data = plot_long,
      aes(x = trend, y = lc_name, fill = scenario),
      shape = 21,
      colour = "grey10",
      stroke = 0.25,
      size = 2.8,
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.5
    ) +
    scale_fill_manual(
      values = c(
        "Unmasked" = "lightblue",
        "Masked" = "darkblue"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = seq(0, trend_xmax, by = 0.001 * scale_factor),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(xlim = c(0, trend_xmax), clip = "off") +
    labs(x = unit_label, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.25),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      plot.caption = element_text(
        size = 7,
        colour = "grey35",
        hjust = 0
      ),
      plot.margin = margin(14, 80, 5.5, 5.5)
    )
  p
}

plot_lc_abs_vs_rel <- function(plot_df) {
  # Generate a scatter plot comparing absolute and relative LAI trends by land
  # cover type using ggplot2. The plot includes points sized by area, colored by
  # vegetation type, and labeled with land cover names. Horizontal and vertical
  # reference lines are added at zero, and the plot is styled with a clean
  # theme.
  veg_colors <- c(
    "Tree Cover" = "#009E73",
    # bluish green
    "Shrubland"  = "#E69F00",
    # orange
    "Herbaceous" = "#56B4E9",
    # sky blue
    "Cropland"   = "#D55E00",
    # vermillion
    "Other"      = "#999999" # grey
  )

  p <- ggplot(plot_df, aes(x = trend_rel, y = trend_abs, color = veg_type)) +
    geom_hline(
      yintercept = 0,
      colour = "grey45",
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = 0,
      colour = "grey45",
      linewidth = 0.3
    ) +
    geom_point(aes(size = area_mkm2), alpha = 0.8, stroke = 0.5) +
    geom_text_repel(
      data = label_df,
      aes(x = trend_rel, y = trend_abs, label = lc_name),
      inherit.aes = FALSE,
      size = 2.7,
      colour = "black",
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 26
    ) +
    scale_color_manual(values = veg_colors, name = "Vegetation Type") +
    scale_size_area(name = "Area (Mkm²)", max_size = 9) +
    labs(
      x = expression("Relative LAI trend (% " * yr^{
        -1
      } * ")"),
      y = expression("Absolute LAI trend (" * m^2 ~ m^{
        -2
      } ~ yr^
        {
          -1
        } * ")"),
      title = "Absolute vs. Relative LAI Trends by Land Cover",
      subtitle = sprintf("Masked trends (%s, %s); Point size = Global area", mask, tau)
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      plot.caption = element_text(size = 8, color = "grey50"),
      plot.margin = margin(10, 15, 10, 10)
    )
  p
}


plot_kg <- function(plot_long, scale_factor) {
  # Generate a plot for kg trends using ggplot2, showing trends and confidence
  # intervals for different kg labels. The plot includes error bars for
  # confidence intervals, points for trends, and a vertical reference line at
  # zero. The x-axis is scaled based on the provided scale_factor.
  finite_x <- c(plot_long$trend, plot_long$ci_lower, plot_long$ci_upper)
  finite_x <- finite_x[is.finite(finite_x)]
  x_min <- 0
  trend_xmax <- max(finite_x) + 0.3 * max(finite_x)
  area_label_x <- max(finite_x) + 0.3 * max(finite_x)
  area_header_x <- area_label_x
  x_range <- trend_xmax - x_min

  p <- ggplot() +
    geom_segment(
      data = plot_tab,
      aes(
        x = mean_unmasked,
        xend = mean_masked,
        y = kg_label,
        yend = kg_label
      ),
      linewidth = 0.7,
      colour = "grey45",
      lineend = "round"
    ) +
    geom_errorbar(
      data = plot_long,
      aes(xmin = ci_lower, xmax = ci_upper, y = kg_label),
      width = 0.12,
      colour = "grey50",
      alpha = 0.5,
      linewidth = 0.4,
      na.rm = TRUE
    ) +
    geom_point(
      data = plot_long,
      aes(
        x = trend,
        y = kg_label,
        fill = scenario,
        shape = factor(shape_type)
      ),
      colour = "grey10",
      stroke = 0.25,
      size = 2.8,
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey30",
      linewidth = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Unmasked" = "lightblue",
        "Masked" = "darkblue"
      ),
      name = NULL
    ) +
    scale_shape_manual(
      values = c(`1` = 1, `16` = 21),
      labels = c("Not significant", "Significant"),
      name = "Significance"
    ) +
    scale_x_continuous(
      breaks = seq(0, trend_xmax, by = 0.001 * scale_factor),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(xlim = c(0, trend_xmax), clip = "off") +
    labs(x = unit_label, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.25),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      plot.caption = element_text(
        size = 7,
        colour = "grey35",
        hjust = 0
      ),
      plot.margin = margin(14, 80, 5.5, 5.5)
    )
}
