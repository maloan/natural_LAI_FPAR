# ==============================================================================
# 09_lai_yearmax_trend_maps.R — Unmasked vs CCI-masked for absolute and relative
# trends
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
source(here("R", "helpers", "plotting.R"))

# config
alpha <- "alpha_0.1"
mask <- "CCI"
var <- "LAI"
metric <- "yearmax"
alpha <- 0.05
include_masked_out <- FALSE

f_p_unm <- here(
  "analysis",
  "unmasked",
  "0p25",
  sprintf("%s_georef_%s_trend_mk_pval_0p25.nc", var, metric)
)
f_p_msk <- here(
  "output",
  alpha,
  "eval",
  sprintf("trend_%s_%s", var, mask),
  sprintf("%s_%s_trend_mk_pval_0p25.nc", var, metric)
)
f_abs_unm <- here(
  "analysis",
  "unmasked",
  "0p25",
  sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, metric)
)
f_abs_msk <- here(
  "output",
  alpha,
  "eval",
  sprintf("trend_%s_%s", var, mask),
  sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, metric)
)
f_rel_unm <- here(
  "analysis",
  "unmasked",
  "0p25",
  sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, metric)
)
f_rel_msk <- here(
  "output",
  alpha,
  "eval",
  sprintf("trend_%s_%s", var, mask),
  sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, metric)
)
outdir_fig <- here("analysis", "results", "figures", "maps")
outdir_tbl <- here("analysis", "results", "tables", "maps")
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tbl, recursive = TRUE, showWarnings = FALSE)
# load rasters
required_files <- c(f_abs_unm, f_abs_msk, f_rel_unm, f_rel_msk, f_p_unm, f_p_msk)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing input raster file(s): ",
    paste(missing_files, collapse = "; ")
  )
}

abs_unm <- rast(f_abs_unm)[[1]]
abs_msk <- rast(f_abs_msk)[[1]]
rel_unm <- rast(f_rel_unm)[[1]]
rel_msk <- rast(f_rel_msk)[[1]]
p_unm <- rast(f_p_unm)[[1]]
p_msk <- rast(f_p_msk)[[1]]
# convert to %/yr:
rel_unm <- 100 * rel_unm
rel_msk <- 100 * rel_msk

# common limits
abs_all <- c(abs_unm, abs_msk)
rel_all <- c(rel_unm, rel_msk)
lims_abs <- sym_q_lim(abs_all, q = 0.99)
lims_abs <- max(abs(lims_abs))
lims_abs <- c(-lims_abs, lims_abs)
lims_rel <- sym_q_lim(rel_all, q = 0.99)
lims_rel <- max(abs(lims_rel))
lims_rel <- c(-lims_rel, lims_rel)

# optional masked-out region
mask_out_abs <- is.na(abs_msk) & !is.na(abs_unm)
mask_out_rel <- is.na(rel_msk) & !is.na(rel_unm)
abs_out <- abs_unm
abs_out[!mask_out_abs] <- NA
rel_out <- rel_unm
rel_out[!mask_out_rel] <- NA

# grey overlay where NOT significant (alpha <- 0.05) but data exists
grey_unm_abs <- (p_unm > alpha) & !is.na(abs_unm)
grey_msk_abs <- (p_msk > alpha) & !is.na(abs_msk)
grey_out_abs <- (p_unm > alpha) & !is.na(abs_out)
grey_unm_rel <- (p_unm > alpha) & !is.na(rel_unm)
grey_msk_rel <- (p_msk > alpha) & !is.na(rel_msk)
grey_out_rel <- (p_unm > alpha) & !is.na(rel_out)
# data frames
df_abs_unm <- to_df(abs_unm, "z")
df_abs_msk <- to_df(abs_msk, "z")
df_abs_out <- to_df(abs_out, "z")
df_rel_unm <- to_df(rel_unm, "z")
df_rel_msk <- to_df(rel_msk, "z")
df_rel_out <- to_df(rel_out, "z")
df_grey_abs_unm <- to_df(grey_unm_abs, "g") |> dplyr::filter(g == 1)
df_grey_abs_msk <- to_df(grey_msk_abs, "g") |> dplyr::filter(g == 1)
df_grey_abs_out <- to_df(grey_out_abs, "g") |> dplyr::filter(g == 1)
df_grey_rel_unm <- to_df(grey_unm_rel, "g") |> dplyr::filter(g == 1)
df_grey_rel_msk <- to_df(grey_msk_rel, "g") |> dplyr::filter(g == 1)
df_grey_rel_out <- to_df(grey_out_rel, "g") |> dplyr::filter(g == 1)

# panel titles
col_titles <- c("Unmasked", "Masked")
fill_abs <- expression("Slope (m"^2 * " m"^{
  -2
} * " yr"^{
  -1
} * ")")
fill_rel <- expression("Trend (% yr"^
  {
    -1
  } * ")")

# build panels
p11 <- plot_map(
  df_abs_unm,
  "z",
  lims_abs,
  title = col_titles[1],
  fill_title = fill_abs,
  df_grey = df_grey_abs_unm
)
p12 <- plot_map(
  df_abs_msk,
  "z",
  lims_abs,
  title = col_titles[2],
  fill_title = fill_abs,
  df_grey = df_grey_abs_msk
)
p13 <- plot_map(
  df_abs_out,
  "z",
  lims_abs,
  title = "Masked-out",
  fill_title = fill_abs,
  df_grey = df_grey_abs_out
)
p21 <- plot_map(
  df_rel_unm,
  "z",
  lims_rel,
  title = col_titles[1],
  fill_title = fill_rel,
  df_grey = df_grey_rel_unm
)
p22 <- plot_map(
  df_rel_msk,
  "z",
  lims_rel,
  title = col_titles[2],
  fill_title = fill_rel,
  df_grey = df_grey_rel_msk
)
p23 <- plot_map(
  df_rel_out,
  "z",
  lims_rel,
  title = "Masked-out",
  fill_title = fill_rel,
  df_grey = df_grey_rel_out
)
if (isTRUE(include_masked_out)) {
  p11 <- add_panel_tag(p11, "(a)")
  p12 <- add_panel_tag(p12, "(b)")
  p13 <- add_panel_tag(p13, "(c)")
  p21 <- add_panel_tag(p21, "(d)")
  p22 <- add_panel_tag(p22, "(e)")
  p23 <- add_panel_tag(p23, "(f)")

  p_abs_grid <- p11 + p12 + p13 + plot_layout(ncol = 3, guides = "collect")
  p_rel_grid <- p21 + p22 + p23 + plot_layout(ncol = 3, guides = "collect")
} else {
  p11 <- add_panel_tag(p11, "(a)")
  p12 <- add_panel_tag(p12, "(b)")
  p21 <- add_panel_tag(p21, "(c)")
  p22 <- add_panel_tag(p22, "(d)")

  p_abs_grid <- p11 + p12 + plot_layout(ncol = 2, guides = "collect")
  p_rel_grid <- p21 + p22 + plot_layout(ncol = 2, guides = "collect")
}

row_abs <- row_label("Absolute slope") + p_abs_grid +
  plot_layout(widths = c(0.06, 1), ncol = 2) &
  theme(legend.position = "bottom")

row_rel <- row_label("Relative trend") + p_rel_grid +
  plot_layout(widths = c(0.06, 1), ncol = 2) &
  theme(legend.position = "bottom")

fig <- row_abs / row_rel

out_png <- file.path(
  outdir_fig,
  sprintf("%s_%s_trend_map_%s_%s_main.png", var, metric, mask, alpha)
)
out_pdf <- file.path(
  outdir_fig,
  sprintf("%s_%s_trend_map_%s_%s_main.pdf", var, metric, mask, alpha)
)
ggsave(out_png,
  fig,
  width = 10.2,
  height = 7.2,
  dpi = 400
)
ggsave(out_pdf, fig, width = 10.2, height = 7.2)
