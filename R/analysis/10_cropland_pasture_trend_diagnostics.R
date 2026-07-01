# ==============================================================================
# 10_cropland_pasture_trend_diagnostics.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(rnaturalearth)
  library(scales)
  library(scico)
  library(here)
})

source(here("R", "helpers", "weighted_means.R"))
source(here("R", "helpers", "plotting.R"))
source(here("R", "helpers", "io.R"))

tau <- "tau_0.1"
var <- "LAI"
metric <- "yearmean"
luh_source <- "CCI"
alpha <- 0.05
limit_q <- 0.95
luh_source <- toupper(luh_source)
if (!luh_source %in% c("CCI", "GLC")) {
  stop("luh_source must be CCI or GLC; got: ", luh_source)
}

tau_num <- as.numeric(sub("^tau_", "", tau))
if (!is.finite(tau_num)) {
  stop("Could not parse tau from config: ", tau)
}
tau_tag2 <- gsub("\\.", "p", sprintf("%.2f", tau_num))
outdir_fig <- here("analysis", "results", "figures", "maps")
outdir_tbl <- here("analysis", "results", "tables", "maps")
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tbl, recursive = TRUE, showWarnings = FALSE)

outcsv <- file.path(outdir_tbl, "cropland_pasture_trend_summary_statistics.csv")
outpng <- file.path(
  outdir_fig,
  "cropland_pasture_trends_absolute_relative_timeseries.png"
)
outpdf <- file.path(
  outdir_fig,
  "cropland_pasture_trends_absolute_relative_timeseries.pdf"
)

f_area005 <- here("src", "area_0p05_validdomain_km2.nc")
f_area025 <- here("src", "area_0p25_validdomain_km2.nc")
f_abs_tr <- here(
  "analysis",
  "unmasked",
  "0p25",
  sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, metric)
)
f_rel_tr <- here(
  "analysis",
  "unmasked",
  "0p25",
  sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, metric)
)
f_p_tr <- here(
  "analysis",
  "unmasked",
  "0p25",
  sprintf("%s_georef_%s_trend_mk_pval_0p25.nc", var, metric)
)
f_abs_tr_cci_nat <- here(
  "output",
  tau,
  "eval",
  sprintf("trend_%s_%s", var, "CCI"),
  sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, metric)
)
f_ts_abs <- here(
  "analysis",
  "unmasked",
  "0p25",
  sprintf("%s_georef_%s_0p25.nc", var, metric)
)
f_mask_cci_005 <- here(
  "output",
  tau,
  "masks",
  "mask_cci",
  sprintf(
    "mask_used_frac_fused_tau%s_k3_1992-2020_0p05.tif",
    tau_tag2
  )
)
f_mask_luh_005 <- here(
  "output",
  tau,
  "masks",
  "mask_luh_overlap",
  sprintf(
    "mask_luh_overlap_%s_Gmin0p10_Pmin0p10_alpha0p50_1992-2020_0p05_rep.tif",
    luh_source
  )
)
required_files <- c(
  f_area005,
  f_area025,
  f_abs_tr,
  f_rel_tr,
  f_p_tr,
  f_abs_tr_cci_nat,
  f_ts_abs,
  f_mask_cci_005,
  f_mask_luh_005
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input file(s): ", paste(missing_files, collapse = "; "))
}

area005 <- rast(f_area005)[[1]]
area025 <- rast(f_area025)[[1]]
abs_tr <- rast(f_abs_tr)[[1]]
rel_tr <- rast(f_rel_tr)[[1]]
p_tr <- rast(f_p_tr)[[1]]
abs_tr_cci_nat <- rast(f_abs_tr_cci_nat)[[1]]
rel_tr <- 100 * rel_tr
ts_abs <- rast(f_ts_abs)
years <- 1982:(1982 + nlyr(ts_abs) - 1)
mask_cci_005 <- rast(f_mask_cci_005)[[1]]
mask_luh_005 <- rast(f_mask_luh_005)[[1]]
mask_frac_to_area025 <- function(mask005, area005) {
  m <- clamp(mask005,
    lower = 0,
    upper = 1,
    values = TRUE
  )
  ok <- is.finite(area005) & area005 > 0 & is.finite(m)
  area_excl005 <- ifel(ok, m * area005, NA)
  aggregate(area_excl005,
    fact = 5,
    fun = "sum",
    na.rm = TRUE
  )
}

mk_map_df <- function(r, w, area025, min_excl_frac = 0) {
  frac_excl <- ifel(is.finite(area025) &
    area025 > 0 & is.finite(w), w / area025, NA)
  z <- ifel(is.finite(frac_excl) &
    frac_excl >= min_excl_frac, r, NA)
  df <- as.data.frame(z, xy = TRUE, na.rm = FALSE)
  names(df) <- c("lon", "lat", "value")
  df
}
mk_nonsig_df <- function(p,
                         w,
                         area025,
                         alpha = 0.05,
                         min_excl_frac = 0) {
  frac_excl <- ifel(is.finite(area025) &
    area025 > 0 & is.finite(w), w / area025, NA)
  z <- ifel(
    is.finite(frac_excl) &
      frac_excl >= min_excl_frac &
      is.finite(p) & p > alpha,
    1,
    NA
  )
  df <- as.data.frame(z, xy = TRUE, na.rm = TRUE)
  names(df) <- c("lon", "lat", "flag")
  df
}

ts_cols <- c(
  "Cropland-excluded (CCI)" = "black",
  "Pasture-overlap (LUH2)" = "black"
)


ts_lm_label <- function(df) {
  d <- df |>
    dplyr::filter(is.finite(.data$year), is.finite(.data$value))

  if (nrow(d) < 3) {
    return("Insufficient data")
  }

  fit <- lm(value ~ year, data = d)
  slope <- unname(coef(fit)[["year"]])
  pval <- summary(fit)$coefficients["year", "Pr(>|t|)"]

  if (abs(slope) < 0.001) {
    slope_str <- sprintf("%.2e", slope)
  } else {
    slope_str <- sprintf("%.4f", slope)
  }

  if (pval < 0.001) {
    p_str <- "p < 0.001"
  } else if (pval < 0.01) {
    p_str <- sprintf("p = %.3f", pval)
  } else {
    p_str <- sprintf("p = %.2f", pval)
  }

  sprintf("Slope: %s %s yr⁻¹\n%s", slope_str, "m2 m-2", p_str)
}

# compute masks on 0.25 deg
w_cci <- mask_frac_to_area025(mask_cci_005, area005)
w_luh <- mask_frac_to_area025(mask_luh_005, area005)

drop_cci_025 <- is.na(abs_tr_cci_nat) & !is.na(abs_tr)
w_cci_map <- ifel(drop_cci_025, area025, NA)

drop_luh_025 <- (is.finite(w_luh) & w_luh > 0) & drop_cci_025
w_luh_map <- ifel(drop_luh_025, area025, NA)

s_cci_abs <- weighted_stats(abs_tr, w_cci_map)
s_cci_rel <- weighted_stats(rel_tr, w_cci_map)
s_luh_abs <- weighted_stats(abs_tr, w_luh_map)
s_luh_rel <- weighted_stats(rel_tr, w_luh_map)

summary_tbl <- tibble(
  region = c("Cropland-excluded (CCI)", "Pasture-overlap (LUH2)"),
  area_km2 = c(s_cci_abs$area, s_luh_abs$area),
  abs_trend_mean = c(s_cci_abs$mean, s_luh_abs$mean),
  abs_trend_sd = c(s_cci_abs$sd, s_luh_abs$sd),
  rel_trend_mean_pct = c(s_cci_rel$mean, s_luh_rel$mean),
  rel_trend_sd_pct = c(s_cci_rel$sd, s_luh_rel$sd)
)
write_csv(round_numeric(summary_tbl, 5), outcsv)

ts_df <- bind_rows(
  tibble(
    year = years,
    value = wmean_series(ts_abs, w_cci_map),
    region = "Cropland-excluded (CCI)"
  ),
  tibble(
    year = years,
    value = wmean_series(ts_abs, w_luh_map),
    region = "Pasture-overlap (LUH2)"
  )
)

ts_ylim <- range(ts_df$value, na.rm = TRUE)
ts_pad <- 0.05 * diff(ts_ylim)
ts_ylim <- ts_ylim + c(-ts_pad, ts_pad)

map_cci_abs <- mk_map_df(abs_tr, w_cci_map, area025, min_excl_frac = 0)
map_cci_rel <- mk_map_df(rel_tr, w_cci_map, area025, min_excl_frac = 0)
map_luh_abs <- mk_map_df(abs_tr, w_luh_map, area025, min_excl_frac = 0)
map_luh_rel <- mk_map_df(rel_tr, w_luh_map, area025, min_excl_frac = 0)
map_cci_nonsig <- mk_nonsig_df(p_tr,
  w_cci_map,
  area025,
  alpha = alpha,
  min_excl_frac = 0
)
map_luh_nonsig <- mk_nonsig_df(p_tr,
  w_luh_map,
  area025,
  alpha = alpha,
  min_excl_frac = 0
)
sym_vec_lim <- function(x, q = 0.95) {
  lim <- as.numeric(stats::quantile(abs(x), probs = q, na.rm = TRUE))
  c(-lim, lim)
}
lims_abs <- sym_vec_lim(c(map_cci_abs$value, map_luh_abs$value), q = limit_q)
lims_rel <- sym_vec_lim(c(map_cci_rel$value, map_luh_rel$value), q = limit_q)
fill_abs <- expression("Slope (" * m^2 ~ m^{
  -2
} ~ yr^
  {
    -1
  } * ")")
fill_rel <- expression("Relative trend (% " * yr^{
  -1
} * ")")

common_legend_theme <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.key.width = unit(2, "cm"),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 7),
  legend.spacing.x = unit(0.1, "cm")
)

p_cci_abs <- plot_map(
  map_cci_abs,
  "value",
  lims_abs,
  title = "CCI-excluded: absolute trend",
  fill_title = fill_abs,
  df_grey = map_cci_nonsig
) +
  guides(fill = guide_colourbar(
    direction = "horizontal",
    title.position = "top",
    barheight = unit(0.25, "cm")
  )) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 0, 1, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, -8, 2)
  )

p_cci_rel <- plot_map(
  map_cci_rel,
  "value",
  lims_rel,
  title = "CCI-excluded: relative trend",
  fill_title = fill_rel,
  df_grey = map_cci_nonsig
) +
  guides(fill = guide_colourbar(
    direction = "horizontal",
    title.position = "top",
    barheight = unit(0.25, "cm")
  )) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 0, 1, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, -8, 2)
  )

p_luh_abs <- plot_map(
  map_luh_abs,
  "value",
  lims_abs,
  title = "LUH overlap: absolute trend",
  fill_title = fill_abs,
  df_grey = map_luh_nonsig
) +
  guides(fill = "none")

p_luh_rel <- plot_map(
  map_luh_rel,
  "value",
  lims_rel,
  title = "LUH overlap: relative trend",
  fill_title = fill_rel,
  df_grey = map_luh_nonsig
) +
  guides(fill = "none")

p_cci_abs <- p_cci_abs + common_legend_theme
p_cci_rel <- p_cci_rel + common_legend_theme
p_luh_abs <- p_luh_abs + common_legend_theme
p_luh_rel <- p_luh_rel + common_legend_theme

fig <- (p_cci_abs + p_cci_rel) / (p_luh_abs + p_luh_rel) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggplot2::ggsave(outpng,
  fig,
  width = 10.2,
  height = 7.2,
  dpi = 400
)
ggplot2::ggsave(outpdf,
  fig,
  width = 10.2,
  height = 7.2,
  device = cairo_pdf
)
