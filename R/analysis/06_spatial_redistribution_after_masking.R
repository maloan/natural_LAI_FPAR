# ==============================================================================
# 06_spatial_redistribution_after_masking.R
# Spatial redistribution after masking
# Map + zonal mean of Δ trend slope: (natural-only − unmasked)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(rnaturalearth)
  library(sf)
  library(scales)
  library(scico)
  library(here)
  library(patchwork)
})

# ---- config ------------------------------------------------------------------
var <- "LAI"
mask <- "CCI"
tau <- "tau_0.1"
metric <- "yearmean"

out_dir <- here(
  "analysis", "results", "spatial_redistribution_after_masking", tau
)

out_dir <- here(
  "analysis", "results", "spatial_redistribution_after_masking", tau
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- theme -------------------------------------------------------------------
theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 13, face = "bold"),
      plot.subtitle    = element_text(size = 10),
      strip.text       = element_text(size = 10, face = "bold"),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 9),
      legend.position  = "left",
      legend.box       = "vertical",
      legend.text      = element_text(size = 9)
    )
}
lat_labels <- function(x) {
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N")))
}
lon_labels <- function(x) {
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°W"), paste0(x, "°E")))
}

lat_breaks_30 <- seq(-90, 90, by = 30)
lon_breaks_60 <- seq(-180, 180, by = 30)


# ---- paths -------------------------------------------------------------------
f_unm <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, metric)
)
f_nat <- here(
  "output", tau, "eval",
  sprintf("trend_%s_%s", var, mask),
  sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, metric)
)
f_area <- here("src", "area_0p25_validdomain_km2.nc")

# ---- load + align ------------------------------------------------------------
b_unm <- rast(f_unm)[[1]]
b_nat <- rast(f_nat)[[1]]
area <- rast(f_area)[[1]]

valid_dom <- is.finite(area) & (area > 0)

# ---- Δ slope raster ----------------------------------------------------------
# Δβ is only defined where both products exist (and within valid domain)
ok_db <- valid_dom & is.finite(b_unm) & is.finite(b_nat)
db <- ifel(ok_db, b_nat - b_unm, NA_real_)

# ---- map limits --------------------------------------------------------------
sym_q_lim <- function(x, q = 0.95, fallback = 1e-3, min_n = 50L) {
  v <- values(x, mat = FALSE)
  v <- v[is.finite(v)]
  if (length(v) < min_n) {
    return(c(-fallback, fallback))
  }
  lim <- as.numeric(stats::quantile(abs(v), probs = q, na.rm = TRUE))
  if (!is.finite(lim) || lim <= 0) {
    lim <- fallback
  }
  c(-lim, lim)
}

lims_db <- sym_q_lim(db)

# ---- zonal mean Δ slope (area-weighted) --------------------------------------
lat <- init(area, "y")
zone <- floor((lat + 90) / 1) + 1

ok_z <- is.finite(db) & is.finite(area) & (area > 0)
num <- ifel(ok_z, db * area, NA)
den <- ifel(ok_z, area, NA)

s_num <- zonal(num, zone, "sum", na.rm = TRUE)
names(s_num) <- c("zone", "num_sum")
s_den <- zonal(den, zone, "sum", na.rm = TRUE)
names(s_den) <- c("zone", "den_sum")

z_db <- merge(s_num, s_den, by = "zone", all = TRUE) |>
  as_tibble() |>
  transmute(
    lat_band = -90 + (zone - 0.5),
    delta_slope = num_sum / den_sum,
    area_sum_km2 = den_sum
  ) |>
  mutate(delta_slope = ifelse(is.finite(delta_slope), delta_slope, NA_real_)) |>
  arrange(lat_band)

out_csv <- file.path(
  out_dir, sprintf("delta_slope_zonal_%s_%s_%s_%s.csv", var, mask, tau, metric)
)
write_csv(z_db, out_csv)

# ---- labels ------------------------------------------------------------------
var_lab <- if (var == "LAI") {
  "LAI"
} else {
  "fAPAR"
}
ylab_map <- sprintf("Δ slope (%s yr\u207B\u00B9)", var_lab)
ylab_zon <- sprintf("Zonal mean Δ slope (%s yr\u207B\u00B9)", var_lab)

title_txt <- sprintf(
  "%s %s: Δ trend slope (natural-only − unmasked)", var_lab, metric
)
subtitle_txt <- sprintf("mask: %s; %s; lat bins = %g°", mask, tau, 1)

# ---- map plot ----------------------------------------------------------------
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

df_db <- as.data.frame(db, xy = TRUE, na.rm = FALSE)
names(df_db) <- c("lon", "lat", "delta_slope")

# Use an explicitly diverging palette for signed Δ fields
div_cols <- scico::scico(256, palette = "vik", direction = 1)
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf") |>
  sf::st_transform(4326)

p_map <- ggplot(df_db) +
  geom_raster(aes(lon, lat, fill = delta_slope), na.rm = FALSE) +
  geom_sf(data = coast, linewidth = 0.2, colour = "black", fill = NA, inherit.aes = FALSE) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  scale_x_continuous(breaks = lon_breaks_60, labels = lon_labels) +
  scale_y_continuous(breaks = lat_breaks_30, labels = lat_labels) +
  scale_fill_gradientn(
    colours = div_cols,
    limits = lims_db,
    oob = squish,
    na.value = "transparent", # makes missing cells not paint white over everything
    name = ylab_map
  ) +
  labs(title = title_txt, subtitle = subtitle_txt, x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank()) +
  theme_pub()


# ---- zonal plot --------------------------------------------------------------
p_zon <- ggplot(z_db, aes(delta_slope, lat_band)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_path(linewidth = 0.85, colour = "black", na.rm = TRUE) +
  scale_y_continuous(breaks = lat_breaks_30, labels = lat_labels) +
  labs(y = NULL, x = ylab_zon) +
  theme_pub() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ---- save outputs ------------------------------------------------------------
out_png_map <- file.path(
  out_dir, sprintf("delta_slope_map_%s_%s_%s_%s.png", var, mask, tau, metric)
)
ggsave(out_png_map, p_map, width = 11, height = 5, dpi = 300)

out_png_zon <- file.path(
  out_dir, sprintf("delta_slope_zonal_%s_%s_%s_%s.png", var, mask, tau, metric)
)
ggsave(out_png_zon, p_zon, width = 4.5, height = 5, dpi = 300)

# combined figure (side-by-side; shared latitude axis for visual comparison)
p_comb <- p_map + p_zon +
  plot_layout(ncol = 2, widths = c(3, 1.0), guides = "collect") &
  theme(legend.position = "left")

out_png_comb <- file.path(
  out_dir, sprintf("delta_slope_map_plus_zonal_%s_%s_%s_%s.png", var, mask, tau, metric)
)
ggsave(out_png_comb, p_comb, width = 14, height = 5.2, dpi = 300)
