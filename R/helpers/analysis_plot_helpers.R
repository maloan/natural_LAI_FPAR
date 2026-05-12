## =============================================================================
## analysis_plot_helpers.R — Shared helpers for analysis plotting
##
## Includes helpers for significance markers, error bars, and theme utilities
## for publication-ready figures with statistical annotations.
## =============================================================================

# Add significance markers to bar/point plots
# Marks CI intervals that don't cross zero with asterisk (*)
geom_sig_marker <- function(data = NULL, mapping = aes(x = x, y = y),
                            sig_col = "sig_flag", y_offset = 0, size = 3.5,
                            color = "black", ...) {
  # Usage: geom_sig_marker(aes(y = value), sig_col = "sig_flag", y_offset = 0.02)
  # Expects data to have sig_flag column with values "" (not sig) or "*" (sig)
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

# Mark non-significant results by fading fill color
geom_sig_col <- function(aes_sig = aes(fill_alpha = sig),
                         alpha_sig = 1.0,
                         alpha_nonsig = 0.3, ...) {
  # Usage: geom_sig_col(aes_sig = aes(fill_alpha = sig_flag == "*"))
  # Fades bars/points where sig_flag != "*"
  scale_alpha_manual(
    values = c("TRUE" = alpha_sig, "FALSE" = alpha_nonsig),
    guide = "none",
    ...
  )
}

# Highlight CIs that cross zero
highlight_ci_crossing <- function(ci_lower, ci_upper,
                                  color_nonsig = "grey70",
                                  color_sig = NA,
                                  linewidth = 1.2) {
  # Returns color vector for error bars based on CI zero-crossing
  # Usage in ggplot: geom_errorbar(aes(color = highlight_ci(...)))
  crosses_zero <- ci_lower * ci_upper <= 0
  ifelse(crosses_zero, color_nonsig, color_sig)
}

#'
#' Computes symmetric limits around zero based on quantile of absolute values.
#' Useful for diverging color scales in trend maps.
#'
#' @param r SpatRaster with values to limit
#' @param q Quantile threshold (0-1), default 0.99
#'
#' @return Numeric vector c(-limit, limit)
#'
sym_q_lim <- function(r, q = 0.99) {
  v <- values(r, mat = FALSE)
  v <- v[is.finite(v)]
  lim <- as.numeric(stats::quantile(abs(v), probs = q, na.rm = TRUE))
  c(-lim, lim)
}

#' Convert raster to data frame with xy coordinates
#'
#' Flattens raster to data frame with lon/lat columns plus value column.
#' Useful for ggplot2 geom_tile mapping.
#'
#' @param r SpatRaster to convert
#' @param name Column name for raster values (default "z")
#'
#' @return Data frame with columns: lon, lat, and name
#'
to_df <- function(r, name = "z") {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df) <- c("lon", "lat", name)
  df
}

#' Format latitude labels with hemisphere
#'
#' Converts latitude values to formatted labels (e.g., "60°N", "30°S", "0°")
#'
#' @param x Numeric latitude values
#' @return Character vector of formatted labels
#'
lat_labels <- function(x) {
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N")))
}

#' Format longitude labels with hemisphere
#'
#' Converts longitude values to formatted labels (e.g., "60°E", "30°W", "0°")
#'
#' @param x Numeric longitude values
#' @return Character vector of formatted labels
#'
lon_labels <- function(x) {
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°W"), paste0(x, "°E")))
}

#' Standard ggplot2 theme for map plots
#'
#' Applies clean map-style theme with cartographic grid and minimal axes.
#' Includes optional legend positioning.
#'
#' @param base_size Base font size (default 10)
#' @param legend_position Position for legend (NULL = no positioning, "bottom", etc.)
#'
#' @return ggplot2 theme object
#'
theme_map <- function(base_size = 10, legend_position = NULL) {
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

#' Generic map plotting function with diverging color scale
#'
#' Creates ggplot2 world map with raster data overlaid, symbolized with
#' diverging color scale. Includes optional coastlines and masking layer.
#'
#' @param df Data frame with lon, lat, and value columns
#' @param zcol Name of value column
#' @param lims Symmetric limits for color scale
#' @param title Plot title
#' @param fill_title Legend title
#' @param df_grey Optional data frame for grey masking overlay
#'
#' @return ggplot object
#'
plot_map <- function(df, zcol, lims, title = NULL, fill_title = NULL, df_grey = NULL) {
  ggplot(df) +
    geom_tile(aes(x = .data$lon, y = .data$lat, fill = .data[[zcol]])) +
    (if (!is.null(df_grey)) {
      geom_tile(
        data = df_grey, inherit.aes = FALSE,
        aes(x = .data$lon, y = .data$lat),
        fill = "grey85", alpha = 1
      )
    }) +
    geom_sf(
      data = coast, # Uses coast from calling environment
      linewidth = 0.15, colour = "black", fill = NA, inherit.aes = FALSE
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-55, 80), expand = FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, by = 60), labels = lon_labels) +
    scale_y_continuous(breaks = c(-60, -30, 0, 30, 60), labels = lat_labels) +
    scale_fill_gradientn(
      colours = scico::scico(256, palette = "bam", direction = 1),
      limits = lims,
      oob = scales::squish,
      na.value = "grey45",
      name = fill_title
    ) +
    labs(title = title) +
    theme_map()
}
