# ==============================================================================
# 04_metrics_mask_effect_zonal.R
# Zonal effect of land-use masking on mean seasonal amplitude (yearamp)
# Computes independent zonal means for each product, then Δ=(natural − unmasked)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(here)
  library(patchwork)
})

# ---- config ------------------------------------------------------------------
var <- "LAI" # "LAI" or "FPAR"
mask <- "CCI" # "CCI" or "GLC"
tau <- "tau_0.1"

band_deg <- 1L # latitude band width in degrees

lab_unm <- "Unmasked (post-abiotic domain)"
lab_nat <- sprintf("Natural-only (%s-based)", mask)
lab_del <- "Δ (natural − unmasked)"

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
      legend.position  = "bottom",
      legend.box       = "vertical",
      legend.text      = element_text(size = 9)
    )
}

# ---- paths -------------------------------------------------------------------
area_path <- here("src", "area_0p25_validdomain_km2.nc")

outdir <- here("analysis", "results", "metrics_mask_effect_zonal", tau)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

unm_path <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_yearamp_0p25.nc", var)
)

msk_path <- here(
  "output", tau, "eval", sprintf("trend_%s_%s", var, mask),
  sprintf("%s_yearamp_0p25.nc", var)
)
area <- rast(area_path)[[1]]

# ---- helpers -----------------------------------------------------------------
lat_labels <- function(x) {
  ifelse(
    x == 0, "0°",
    ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N"))
  )
}

# time-mean area raster that may have multiple layers (annual time series)
time_mean <- function(r) {
  if (nlyr(r) == 1) {
    return(r[[1]])
  }
  app(r, \(x) mean(x, na.rm = TRUE))
}

# independent zonal area-weighted mean by latitude bands
zonal_wmean_latbands <- function(r, area, band_deg = 1L) {
  compareGeom(area, r, stopOnError = TRUE)
  lat <- init(r, "y")
  zone <- floor((lat + 90) / band_deg) + 1
  ok <- is.finite(r) & is.finite(area) & (area > 0)

  num <- ifel(ok, r * area, NA)
  den <- ifel(ok, area, NA)

  s_num <- zonal(num, zone, "sum", na.rm = TRUE)
  names(s_num) <- c("zone", "num")
  s_den <- zonal(den, zone, "sum", na.rm = TRUE)
  names(s_den) <- c("zone", "den")

  merge(s_num, s_den, by = "zone", all = TRUE) |>
    as_tibble() |>
    transmute(
      lat_band = -90 + (zone - 0.5) * band_deg,
      mean_yearamp = round((num / den), 4),
      area_km2 = den
    ) |>
    arrange(.data$lat_band)
}

# ---- load + compute ----------------------------------------------------------
r_unm <- rast(unm_path)
r_nat <- rast(msk_path)

u_amp <- time_mean(r_unm)
n_amp <- time_mean(r_nat)

z_unm <- zonal_wmean_latbands(u_amp, area, band_deg = band_deg) |>
  mutate(domain = lab_unm)

z_nat <- zonal_wmean_latbands(n_amp, area, band_deg = band_deg) |>
  mutate(domain = lab_nat)

# Δ of zonal means
z_del <- z_unm |>
  select(lat_band, mean_yearamp_unm = mean_yearamp, area_unm_km2 = area_km2) |>
  left_join(
    z_nat |>
      select(
        lat_band,
        mean_yearamp_nat = mean_yearamp, area_nat_km2 = area_km2
      ),
    by = "lat_band"
  ) |>
  transmute(
    lat_band,
    domain = lab_del,
    mean_yearamp = mean_yearamp_nat - mean_yearamp_unm,
    area_unm_km2,
    area_nat_km2
  )

zonal_tbl <- bind_rows(
  z_unm |> transmute(lat_band, domain, mean_yearamp, area_km2,
    area_unm_km2 = area_km2, area_nat_km2 = NA_real_
  ),
  z_nat |> transmute(lat_band, domain, mean_yearamp, area_km2,
    area_unm_km2 = NA_real_, area_nat_km2 = area_km2
  ),
  z_del |> transmute(lat_band, domain, mean_yearamp,
    area_km2 = NA_real_, area_unm_km2, area_nat_km2
  )
)
# ---- write table -------------------------------------------------------------
out_csv <- file.path(
  outdir,
  sprintf("zonal_yearamp_timeMean_%s_%s_%s.csv", var, mask, tau)
)
write_csv(mutate(zonal_tbl, across(where(is.double), ~ round(.x, 4))), out_csv)

# contributing areas for Δ
out_csv2 <- file.path(
  outdir,
  sprintf("zonal_yearamp_timeMean_deltaAreas_%s_%s_%s.csv", var, mask, tau)
)
write_csv(mutate(z_del, across(where(is.double), ~ round(.x, 4))), out_csv2)

# ---- plot --------------------------------------------------------------------
z_abs <- zonal_tbl |>
  filter(domain %in% c(lab_unm, lab_nat), is.finite(mean_yearamp))
z_d <- zonal_tbl |>
  filter(domain == lab_del, is.finite(mean_yearamp))

col_abs <- c(setNames("grey20", lab_unm), setNames("#2C7FB8", lab_nat))

ylab_abs <- if (toupper(var) == "LAI") {
  "Mean seasonal amplitude (LAI)"
} else {
  "Mean seasonal amplitude (fAPAR)"
}
ylab_del <- if (toupper(var) == "LAI") {
  "Δ seasonal amplitude (LAI)"
} else {
  "Δ seasonal amplitude (fAPAR)"
}

p_abs <- ggplot(z_abs, aes(lat_band, mean_yearamp, colour = domain)) +
  geom_line(linewidth = 0.85, na.rm = TRUE) +
  scale_x_continuous(
    limits = c(-90, 90),
    breaks = seq(-90, 90, by = 10),
    labels = lat_labels
  ) +
  scale_colour_manual(values = col_abs) +
  labs(x = NULL, y = ylab_abs, colour = NULL) +
  theme_pub()

p_del <- ggplot(z_d, aes(lat_band, mean_yearamp)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_area(fill = "#D55E00", alpha = 0.25, na.rm = TRUE) +
  geom_line(colour = "#D55E00", linewidth = 0.85, na.rm = TRUE) +
  scale_x_continuous(
    limits = c(-90, 90),
    breaks = seq(-90, 90, by = 10),
    labels = lat_labels
  ) +
  labs(x = "Latitude", y = ylab_del) +
  theme_pub() +
  theme(legend.position = "none")

p <- p_abs / p_del + plot_layout(heights = c(1.15, 1.0)) +
  plot_annotation(
    title = sprintf(
      "%s: mean seasonal amplitude before/after masking (%s; %s)",
      var, mask, tau
    ),
    subtitle = sprintf(
      "Zonal means computed independently per product; Δ is the difference of zonal means (band = %d°).",
      band_deg
    )
  )

out_png <- file.path(
  outdir,
  sprintf("zonal_yearamp_timeMean_%s_%s_%s.png", var, mask, tau)
)
ggsave(out_png, p, width = 11, height = 7.2, dpi = 300)
