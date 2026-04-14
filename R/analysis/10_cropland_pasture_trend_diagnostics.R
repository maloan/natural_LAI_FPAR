# ==============================================================================
# 10_cropland_pasture_trend_diagnostics.R
# Figure 5: Trends in excluded cropland and pasture-influenced regions
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

# ---- config ------------------------------------------------------------------
tau <- "tau_0.1"
var <- "LAI"
metric <- "yearmean"
alpha <- 0.1
limit_q <- 0.95

tau_num <- as.numeric(sub("^tau_", "", tau))
if (!is.finite(tau_num)) {
  stop("Could not parse tau from config: ", tau)
}
tau_tag2 <- gsub("\\.", "p", sprintf("%.2f", tau_num))

outdir <- here("analysis", "results", "paper_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outcsv <- file.path(outdir, "cropland_pasture_trend_summary_statistics.csv")
outpng <- file.path(outdir, "cropland_pasture_trends_absolute_relative_timeseries.png")
outpdf <- file.path(outdir, "cropland_pasture_trends_absolute_relative_timeseries.pdf")

# ---- inputs ------------------------------------------------------------------
f_area005 <- here("src", "area_0p05_validdomain_km2.nc")
f_area025 <- here("src", "area_0p25_validdomain_km2.nc")
f_abs_tr <- here("analysis", "unmasked", "0p25", sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, metric))
f_rel_tr <- here("analysis", "unmasked", "0p25", sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, metric))
f_p_tr <- here("analysis", "unmasked", "0p25", sprintf("%s_georef_%s_trend_mk_pval_0p25.nc", var, metric))
f_abs_tr_cci_nat <- here(
  "output", tau, "eval", sprintf("trend_%s_%s", var, "CCI"),
  sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, metric)
)

f_ts_abs <- here("analysis", "unmasked", "0p25", sprintf("%s_georef_%s_0p25.nc", var, metric))
f_mask_cci_005 <- here(
  "output", tau, "masks", "mask_cci",
  sprintf("mask_used_frac_fused_tau%s_k3_1992-2020_0p05.tif", tau_tag2)
)
f_mask_luh_005 <- here(
  "output", tau, "masks", "mask_luh_overlap",
  sprintf("mask_luh_overlap_CCI_Gmin%s_Pmin%s_alpha0p50_1992-2020_0p05_rep.tif", tau_tag2, tau_tag2)
)

required_files <- c(
  f_area005, f_area025, f_abs_tr, f_rel_tr, f_p_tr, f_abs_tr_cci_nat,
  f_ts_abs, f_mask_cci_005, f_mask_luh_005
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

compareGeom(mask_cci_005, area005, stopOnError = TRUE)
compareGeom(mask_luh_005, area005, stopOnError = TRUE)
compareGeom(abs_tr, area025, stopOnError = TRUE)
compareGeom(rel_tr, area025, stopOnError = TRUE)
compareGeom(p_tr, area025, stopOnError = TRUE)
compareGeom(abs_tr_cci_nat, area025, stopOnError = TRUE)
compareGeom(ts_abs, area025, stopOnError = TRUE)

# ---- helpers -----------------------------------------------------------------
lat_labels <- function(x) ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N")))
lon_labels <- function(x) ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°W"), paste0(x, "°E")))

lon_breaks_60 <- seq(-180, 180, by = 60)

# Subset breaks for cropped map (55S to 80N)
lat_breaks_crop <- c(-60, -30, 0, 30, 60)

ts_cols <- setNames(
  scico::scico(2, palette = "bam", direction = 1),
  c("Cropland-excluded (CCI)", "Pasture-overlap (LUH2)")
)

div_cols <- scico::scico(256, palette = "bam", direction = 1)
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf") |>
  sf::st_transform(4326)

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
      legend.key.height = unit(0.5, "cm"),
      legend.position = "bottom"
    )
}

mask_frac_to_area025 <- function(mask005, area005) {
  m <- clamp(mask005, lower = 0, upper = 1, values = TRUE)
  ok <- is.finite(area005) & area005 > 0 & is.finite(m)
  area_excl005 <- ifel(ok, m * area005, NA)
  aggregate(area_excl005, fact = 5, fun = "sum", na.rm = TRUE)
}

wmean_series <- function(r, w) {
  out <- numeric(nlyr(r))
  for (i in seq_len(nlyr(r))) {
    x <- r[[i]]
    ok <- is.finite(x) & is.finite(w) & w > 0
    num <- as.numeric(global(ifel(ok, x * w, NA), "sum", na.rm = TRUE)[1, 1])
    den <- as.numeric(global(ifel(ok, w, NA), "sum", na.rm = TRUE)[1, 1])
    out[i] <- ifelse(is.finite(den) && den > 0, num / den, NA_real_)
  }
  out
}

weighted_stats <- function(x, w) {
  dx <- as.data.frame(c(x, w), na.rm = FALSE)
  names(dx) <- c("x", "w")
  dx <- dx |> filter(is.finite(x), is.finite(w), w > 0)
  if (nrow(dx) == 0) {
    return(list(mean = NA_real_, sd = NA_real_, area = NA_real_))
  }
  mu <- weighted.mean(dx$x, dx$w)
  vv <- sum(dx$w * (dx$x - mu)^2) / sum(dx$w)
  list(mean = mu, sd = sqrt(vv), area = sum(dx$w))
}

mk_map_df <- function(r, w, area025, min_excl_frac = 0) {
  frac_excl <- ifel(is.finite(area025) & area025 > 0 & is.finite(w), w / area025, NA)
  z <- ifel(is.finite(frac_excl) & frac_excl >= min_excl_frac, r, NA)
  df <- as.data.frame(z, xy = TRUE, na.rm = FALSE)
  names(df) <- c("lon", "lat", "value")
  df
}

mk_nonsig_df <- function(p, w, area025, alpha = 0.1, min_excl_frac = 0) {
  frac_excl <- ifel(is.finite(area025) & area025 > 0 & is.finite(w), w / area025, NA)
  z <- ifel(is.finite(frac_excl) & frac_excl >= min_excl_frac & is.finite(p) & p > alpha, 1, NA)
  df <- as.data.frame(z, xy = TRUE, na.rm = TRUE)
  names(df) <- c("lon", "lat", "flag")
  df
}

sym_lims <- function(x, q = 0.95) {
  lim <- as.numeric(stats::quantile(abs(x), probs = q, na.rm = TRUE))
  c(-lim, lim)
}

plot_map <- function(df, zcol, lims, title = NULL, fill_title = NULL, df_grey = NULL) {
  ggplot(df) +
    geom_tile(aes(x = .data$lon, y = .data$lat, fill = .data[[zcol]])) +
    (if (!is.null(df_grey)) {
      geom_tile(
        data = df_grey, inherit.aes = FALSE,
        aes(x = .data$lon, y = .data$lat),
        fill = "grey85", alpha = 1  # Light grey for non-significant
      )
    }) +
    geom_sf(
      data = coast,
      linewidth = 0.15, colour = "black", fill = NA, inherit.aes = FALSE
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-55, 80), expand = FALSE) +
    scale_x_continuous(breaks = lon_breaks_60, labels = lon_labels) +
    scale_y_continuous(breaks = lat_breaks_crop, labels = lat_labels) +
    scale_fill_gradientn(
      colours = div_cols,
      limits = lims,
      oob = squish,
      na.value = "grey45",  # Medium grey for masked regions
      name = fill_title
    ) +
    labs(title = title) +
    theme_map()
}

theme_ts <- function() {
  theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      plot.title = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.position = "bottom"
    )
}

plot_ts <- function(df, ttl) {
  ggplot(df, aes(.data$year, .data$value, colour = .data$region)) +
    geom_line(linewidth = 0.75, na.rm = TRUE) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, linetype = "dashed", na.rm = TRUE) +
    scale_colour_manual(values = ts_cols) +
    scale_x_continuous(breaks = seq(1982, 2024, by = 7)) +
    labs(title = ttl, x = "Year", y = "Global mean absolute LAI", colour = NULL) +
    theme_ts()
}

# ---- compute masks on 0.25 deg ------------------------------------------------
w_cci <- mask_frac_to_area025(mask_cci_005, area005)
w_luh <- mask_frac_to_area025(mask_luh_005, area005)

# Use exact CCI dropped footprint (vs natural CCI map) for plotting to avoid overlap.
drop_cci_025 <- is.na(abs_tr_cci_nat) & !is.na(abs_tr)
w_cci_map <- ifel(drop_cci_025, area025, NA)

# Use LUH overlap footprint constrained to CCI-dropped cells to avoid overlap with natural CCI area.
drop_luh_025 <- (is.finite(w_luh) & w_luh > 0) & drop_cci_025
w_luh_map <- ifel(drop_luh_025, area025, NA)

compareGeom(w_cci, area025, stopOnError = TRUE)
compareGeom(w_luh, area025, stopOnError = TRUE)
compareGeom(w_cci_map, area025, stopOnError = TRUE)
compareGeom(w_luh_map, area025, stopOnError = TRUE)

# ---- summary stats ------------------------------------------------------------
s_cci_abs <- weighted_stats(abs_tr, w_cci)
s_cci_rel <- weighted_stats(rel_tr, w_cci)
s_luh_abs <- weighted_stats(abs_tr, w_luh)
s_luh_rel <- weighted_stats(rel_tr, w_luh)

summary_tbl <- tibble(
  region = c("Cropland-excluded (CCI)", "Pasture-overlap (LUH2)"),
  area_km2 = c(s_cci_abs$area, s_luh_abs$area),
  abs_trend_mean = c(s_cci_abs$mean, s_luh_abs$mean),
  abs_trend_sd = c(s_cci_abs$sd, s_luh_abs$sd),
  rel_trend_mean_pct = c(s_cci_rel$mean, s_luh_rel$mean),
  rel_trend_sd_pct = c(s_cci_rel$sd, s_luh_rel$sd)
)
write_csv(summary_tbl |> mutate(across(where(is.double), ~ round(.x, 5))), outcsv)

# ---- timeseries ---------------------------------------------------------------
ts_df <- bind_rows(
  tibble(year = years, value = wmean_series(ts_abs, w_cci), region = "Cropland-excluded (CCI)"),
  tibble(year = years, value = wmean_series(ts_abs, w_luh), region = "Pasture-overlap (LUH2)")
)

# ---- maps ---------------------------------------------------------------------
map_cci_abs <- mk_map_df(abs_tr, w_cci_map, area025, min_excl_frac = 0)
map_cci_rel <- mk_map_df(rel_tr, w_cci_map, area025, min_excl_frac = 0)
map_luh_abs <- mk_map_df(abs_tr, w_luh_map, area025, min_excl_frac = 0)
map_luh_rel <- mk_map_df(rel_tr, w_luh_map, area025, min_excl_frac = 0)

map_cci_nonsig <- mk_nonsig_df(p_tr, w_cci_map, area025, alpha = alpha, min_excl_frac = 0)
map_luh_nonsig <- mk_nonsig_df(p_tr, w_luh_map, area025, alpha = alpha, min_excl_frac = 0)

lims_abs <- sym_lims(c(map_cci_abs$value, map_luh_abs$value), q = limit_q)
lims_rel <- sym_lims(c(map_cci_rel$value, map_luh_rel$value), q = limit_q)

fill_abs <- "Slope (LAI yr^-1)"
fill_rel <- "Trend (% yr^-1)"

p1 <- plot_map(map_cci_abs, "value", lims_abs, title = "Cropland-excluded: absolute trend", fill_title = fill_abs, df_grey = map_cci_nonsig)
p2 <- plot_map(map_cci_rel, "value", lims_rel, title = "Cropland-excluded: relative trend", fill_title = fill_rel, df_grey = map_cci_nonsig)
p3 <- plot_ts(ts_df |> filter(region == "Cropland-excluded (CCI)"), "Cropland-excluded: absolute LAI time series")

p4 <- plot_map(map_luh_abs, "value", lims_abs, title = "Pasture-overlap: absolute trend", fill_title = fill_abs, df_grey = map_luh_nonsig)
p5 <- plot_map(map_luh_rel, "value", lims_rel, title = "Pasture-overlap: relative trend", fill_title = fill_rel, df_grey = map_luh_nonsig)
p6 <- plot_ts(ts_df |> filter(region == "Pasture-overlap (LUH2)"), "Pasture-overlap: absolute LAI time series")

fig <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(
    title = "Greening signal within excluded cropland and pasture-influenced regions",
    subtitle = "Top row: CCI cropland-excluded footprint (disjoint from natural CCI); bottom row: LUH2 pasture-overlap footprint constrained to the CCI-dropped area"
  )

ggsave(outpng, fig, width = 13.2, height = 10.2, dpi = 400)
ggsave(outpdf, fig, width = 13.2, height = 10.2)
