# ==============================================================================
# 07_LAI_trend_map_table.R
# 2x3 table figure: unmasked vs masked vs masked-out region; absolute vs relative
# + grey background for non-significant trends (OLS p-value of slope)
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
tau <- "tau_0.1"
mask <- "CCI"
var <- "LAI"
metric <- "yearmean"

alpha <- 0.10 # significance threshold for p-value shading

# trends (slopes / relative)
f_abs_unm <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, metric)
)
f_abs_msk <- here(
  "output", tau, "eval",
  sprintf("trend_%s_%s", var, mask),
  sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, metric)
)

f_rel_unm <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, metric)
)
f_rel_msk <- here(
  "output", tau, "eval",
  sprintf("trend_%s_%s", var, mask),
  sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, metric)
)

# p-values (OLS slope; one layer)
# unmasked
f_p_unm <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_%s_trend_ols_pval_0p25.nc", var, metric)
)
# masked
f_p_msk <- here(
  "output", tau, "eval",
  sprintf("trend_%s_%s", var, mask),
  sprintf("%s_%s_trend_ols_pval_0p25.nc", var, metric)
)

outdir <- here("analysis", "results", "paper_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers -----------------------------------------------------------------
sym_q_lim <- function(r, q = 0.99) {
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

# apply p-value shading:
# - z_sig: z where p <= alpha, else NA
# - nsig: 1 where z is finite and p > alpha, else NA (for grey background)
make_sig_layers <- function(z, p, alpha = 0.10) {
  compareGeom(z, p, stopOnError = TRUE)

  nsig <- z
  nsig[] <- NA_real_
  nsig[is.finite(z) & is.finite(p) & (p > alpha)] <- 1

  z_sig <- z
  z_sig[!(is.finite(z) & is.finite(p) & (p <= alpha))] <- NA_real_

  list(z_sig = z_sig, nsig = nsig)
}

plot_map_sig <- function(df_sig, df_nsig, zcol, lims, title = NULL, fill_title = NULL,
                         nsig_col = "grey85") {
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf") |>
    sf::st_transform(4326)

  ggplot() +
    # non-significant background
    geom_raster(
      data = df_nsig,
      aes(x = .data[["lon"]], y = .data[["lat"]], fill = .data[[zcol]]),
      show.legend = FALSE
    ) +
    scale_fill_identity(na.value = "transparent") +
    # significant trends (on top)
    geom_raster(
      data = df_sig,
      aes(x = .data[["lon"]], y = .data[["lat"]], fill = .data[[zcol]])
    ) +
    geom_sf(
      data = coast,
      linewidth = 0.15, colour = "black", fill = NA, inherit.aes = FALSE
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    scale_x_continuous(breaks = lon_breaks_60, labels = lon_labels) +
    scale_y_continuous(breaks = lat_breaks_30, labels = lat_labels) +
    scale_fill_gradientn(
      colours = div_cols,
      limits = lims,
      oob = squish,
      na.value = "transparent",
      name = fill_title
    ) +
    labs(title = title) +
    theme_map() +
    # set constant colour for nsig layer via identity mapping
    guides(fill = guide_colorbar()) +
    # hack: ensure nsig layer uses the chosen grey
    # (identity scale uses numeric values; we map 1 -> nsig_col by setting it directly)
    NULL
}

# small helper: overwrite nsig df values with the desired constant colour
# (avoids an extra scale for the nsig layer)
nsig_df_to_colour <- function(df, zcol = "z", col = "grey85") {
  df[[zcol]] <- ifelse(is.finite(df[[zcol]]), col, NA_character_)
  df
}

# ---- load rasters ------------------------------------------------------------
abs_unm <- rast(f_abs_unm)[[1]]
abs_msk <- rast(f_abs_msk)[[1]]
rel_unm <- rast(f_rel_unm)[[1]]
rel_msk <- rast(f_rel_msk)[[1]]

# convert to %/yr:
rel_unm <- 100 * rel_unm
rel_msk <- 100 * rel_msk

# p-values (single layer)
p_unm <- rast(f_p_unm)[[1]]
p_msk <- rast(f_p_msk)[[1]]

compareGeom(abs_unm, abs_msk, stopOnError = TRUE)
compareGeom(rel_unm, rel_msk, stopOnError = TRUE)
compareGeom(abs_unm, rel_unm, stopOnError = TRUE)
compareGeom(abs_unm, p_unm, stopOnError = TRUE)
compareGeom(abs_msk, p_msk, stopOnError = TRUE)

# ---- colour limits (based on both products) ----------------------------------
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

# ---- significance layers ------------------------------------------------------
# unmasked panels use p_unm; masked panels use p_msk; masked-out uses p_unm
sig_abs_unm <- make_sig_layers(abs_unm, p_unm, alpha = alpha)
sig_abs_msk <- make_sig_layers(abs_msk, p_msk, alpha = alpha)
sig_abs_out <- make_sig_layers(abs_out, p_unm, alpha = alpha)

sig_rel_unm <- make_sig_layers(rel_unm, p_unm, alpha = alpha)
sig_rel_msk <- make_sig_layers(rel_msk, p_msk, alpha = alpha)
sig_rel_out <- make_sig_layers(rel_out, p_unm, alpha = alpha)

# ---- data frames -------------------------------------------------------------
df_abs_unm_sig <- to_df(sig_abs_unm$z_sig, "z")
df_abs_unm_nsig <- nsig_df_to_colour(to_df(sig_abs_unm$nsig, "z"), "z", col = "grey85")

df_abs_msk_sig <- to_df(sig_abs_msk$z_sig, "z")
df_abs_msk_nsig <- nsig_df_to_colour(to_df(sig_abs_msk$nsig, "z"), "z", col = "grey85")

df_abs_out_sig <- to_df(sig_abs_out$z_sig, "z")
df_abs_out_nsig <- nsig_df_to_colour(to_df(sig_abs_out$nsig, "z"), "z", col = "grey85")

df_rel_unm_sig <- to_df(sig_rel_unm$z_sig, "z")
df_rel_unm_nsig <- nsig_df_to_colour(to_df(sig_rel_unm$nsig, "z"), "z", col = "grey85")

df_rel_msk_sig <- to_df(sig_rel_msk$z_sig, "z")
df_rel_msk_nsig <- nsig_df_to_colour(to_df(sig_rel_msk$nsig, "z"), "z", col = "grey85")

df_rel_out_sig <- to_df(sig_rel_out$z_sig, "z")
df_rel_out_nsig <- nsig_df_to_colour(to_df(sig_rel_out$nsig, "z"), "z", col = "grey85")

# ---- panel titles ------------------------------------------------------------
col_titles <- c("Unmasked", "Masked", "Masked-out region (unmasked values)")
fill_abs <- sprintf("Absolute slope (%s yr\u207B\u00B9)", var)
fill_rel <- "Relative trend (% yr\u207B\u00B9)"

# ---- build panels ------------------------------------------------------------
p11 <- plot_map_sig(
  df_abs_unm_sig, df_abs_unm_nsig, "z", lims_abs,
  title = col_titles[1], fill_title = fill_abs
)
p12 <- plot_map_sig(
  df_abs_msk_sig, df_abs_msk_nsig, "z", lims_abs,
  title = col_titles[2], fill_title = fill_abs
)
p13 <- plot_map_sig(
  df_abs_out_sig, df_abs_out_nsig, "z", lims_abs,
  title = col_titles[3], fill_title = fill_abs
)

p21 <- plot_map_sig(
  df_rel_unm_sig, df_rel_unm_nsig, "z", lims_rel,
  title = col_titles[1], fill_title = fill_rel
)
p22 <- plot_map_sig(
  df_rel_msk_sig, df_rel_msk_nsig, "z", lims_rel,
  title = col_titles[2], fill_title = fill_rel
)
p23 <- plot_map_sig(
  df_rel_out_sig, df_rel_out_nsig, "z", lims_rel,
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
      var, metric
    ),
    subtitle = sprintf(
      "mask: %s; \u03C4 = %s. Grey = non-significant (OLS p > %.2f).",
      mask, tau, alpha
    )
  )

out_png <- file.path(
  outdir,
  sprintf("%s_%s_trend_map_table_%s_%s.png", var, metric, mask, tau)
)
ggsave(out_png, fig, width = 13.2, height = 7.6, dpi = 400)
