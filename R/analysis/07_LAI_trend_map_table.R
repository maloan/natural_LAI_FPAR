# ==============================================================================
# 07_LAI_trend_map_table.R
# 2x3 table figure: unmasked vs masked vs delta; absolute vs relative
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(scales)
  library(scico)
  library(patchwork)
  library(here)
  library(grid)
})

# ---- config ------------------------------------------------------------------
TAU <- "tau_0.1"
MASK <- "CCI"
VAR <- "LAI"
METRIC <- "yearmean"

ROOT <- here::here()

f_abs_unm <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", VAR, METRIC)
)
f_abs_msk <- here(
  "output", TAU, "eval",
  sprintf("trend_%s_%s", VAR, MASK),
  sprintf("%s_%s_trend_slope_peryear_0p25.nc", VAR, METRIC)
)

f_rel_unm <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", VAR, METRIC)
)
f_rel_msk <- here(
  "output", TAU, "eval",
  sprintf("trend_%s_%s", VAR, MASK),
  sprintf("%s_%s_trend_relative_peryear_0p25.nc", VAR, METRIC)
)

OUTDIR <- here("analysis", "results", "paper_figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---- helpers -----------------------------------------------------------------
sym_q_lim <- function(r, q = 0.99, min_n = 50L) {
  v <- values(r, mat = FALSE)
  v <- v[is.finite(v)]
  lim <- as.numeric(stats::quantile(abs(v), probs = q, na.rm = TRUE))
  c(-lim, lim)
}

to_df <- function(r, name = "z") {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df) <- c("lon", "lat", name)
  df
}

lat_labels <- function(x) {
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N")))
}
lon_labels <- function(x) {
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°W"), paste0(x, "°E")))
}
lat_breaks_30 <- seq(-90, 90, by = 30)
lon_breaks_60 <- seq(-180, 180, by = 60)

div_cols <- scico::scico(256, palette = "bam", direction = 1)

theme_map <- function() {
  theme_bw(base_size = 10) +
    theme(
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
}

plot_map <- function(df, zcol, lims, title = NULL, fill_title = NULL) {
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf") |>
    sf::st_transform(4326)

  ggplot(df) +
    geom_raster(aes(x = lon, y = lat, fill = .data[[zcol]])) +
    geom_sf(
      data = coast,
      linewidth = 0.15, colour = "black", fill = NA, inherit.aes = FALSE
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    scale_x_continuous(
      breaks = lon_breaks_60,
      labels = lon_labels
    ) +
    scale_y_continuous(
      breaks = lat_breaks_30,
      labels = lat_labels
    ) +
    scale_fill_gradientn(
      colours = div_cols,
      limits = lims,
      oob = squish,
      na.value = "transparent",
      name = fill_title
    ) +
    labs(title = title) +
    theme_map()
}

# ---- load rasters ------------------------------------------------------------
abs_unm <- rast(f_abs_unm)[[1]]
abs_msk <- rast(f_abs_msk)[[1]]
rel_unm <- rast(f_rel_unm)[[1]]
rel_msk <- rast(f_rel_msk)[[1]]

# convert to %/yr:
rel_unm <- 100 * rel_unm
rel_msk <- 100 * rel_msk

compareGeom(abs_unm, abs_msk, stopOnError = TRUE)
compareGeom(rel_unm, rel_msk, stopOnError = TRUE)
compareGeom(abs_unm, rel_unm, stopOnError = TRUE)

# ---- deltas ------------------------------------------------------------------
abs_all <- c(abs_unm, abs_msk)
rel_all <- c(rel_unm, rel_msk)

lims_abs <- sym_q_lim(abs_all, q = 0.95)
lims_abs <- max(abs(lims_abs))
lims_abs <- c(-lims_abs, lims_abs)

lims_rel <- sym_q_lim(rel_all, q = 0.95)
lims_rel <- max(abs(lims_rel))
lims_rel <- c(-lims_rel, lims_rel)

# ---- masked-out region (where masked is NA but unmasked is not) --------------
mask_out_abs <- is.na(abs_msk) & !is.na(abs_unm)
mask_out_rel <- is.na(rel_msk) & !is.na(rel_unm)

abs_out <- abs_unm
abs_out[!mask_out_abs] <- NA

rel_out <- rel_unm
rel_out[!mask_out_rel] <- NA

# ---- data frames -------------------------------------------------------------
df_abs_unm <- to_df(abs_unm, "z")
df_abs_msk <- to_df(abs_msk, "z")
df_abs_out <- to_df(abs_out, "z")

df_rel_unm <- to_df(rel_unm, "z")
df_rel_msk <- to_df(rel_msk, "z")
df_rel_out <- to_df(rel_out, "z")

# ---- panel titles ------------------------------------------------------------
col_titles <- c("Unmasked", "Masked", "Masked-out region (unmasked values)")
fill_abs <- sprintf("Absolute slope (%s yr\u207B\u00B9)", VAR)
fill_rel <- "Relative trend (% yr\u207B\u00B9)"


# ---- build panels ------------------------------------------------------------
p11 <- plot_map(
  df_abs_unm, "z", lims_abs,
  title = col_titles[1], fill_title = fill_abs
)
p12 <- plot_map(
  df_abs_msk, "z", lims_abs,
  title = col_titles[2], fill_title = fill_abs
)
p13 <- plot_map(
  df_abs_out, "z", lims_abs,
  title = col_titles[3], fill_title = fill_abs
)

p21 <- plot_map(
  df_rel_unm, "z", lims_rel,
  title = col_titles[1], fill_title = fill_rel
)
p22 <- plot_map(
  df_rel_msk, "z", lims_rel,
  title = col_titles[2], fill_title = fill_rel
)
p23 <- plot_map(
  df_rel_out, "z", lims_rel,
  title = col_titles[3], fill_title = fill_rel
)


row_label <- function(txt) {
  wrap_elements(
    full = textGrob(txt, rot = 90, gp = gpar(fontface = "bold", fontsize = 11))
  )
}

row_abs <- row_label("Absolute slope") + (p11 + p12 + p13) +
  plot_layout(widths = c(0.06, 1), ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

row_rel <- row_label("Relative trend") + (p21 + p22 + p23) +
  plot_layout(widths = c(0.06, 1), ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

fig <- (row_abs / row_rel) +
  plot_annotation(
    title = sprintf(
      "%s %s: trend maps (unmasked vs masked and masked-out region)",
      VAR, METRIC
    ),
    subtitle = sprintf("Mask: %s; τ = %s.", MASK, TAU)
  )

out_png <- file.path(
  OUTDIR,
  sprintf("%s_%s_trend_map_table_%s_%s.png", VAR, METRIC, MASK, TAU)
)
ggsave(out_png, fig, width = 13.2, height = 7.6, dpi = 400)
