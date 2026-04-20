# ==============================================================================
# 05_zonal_seasonal_amplitude.R
# Figure 4 component: zonal mean seasonal amplitude (absolute values)
# All scenarios overlaid (unmasked, CCI tau sweep, GLC)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(here)
})

source(here("R", "helpers", "weighted_means.R"))
source(here("R", "helpers", "plotting.R"))

# ---- config ------------------------------------------------------------------
var <- "LAI" # "LAI" or "FPAR"
taus_cci <- c("tau_0.05", "tau_0.1", "tau_0.2")
tau_glc <- "tau_0.1"
tau_ref <- "tau_0.1"

band_deg <- 1L # latitude band width in degrees

scenario_levels <- c("Unmasked", "CCI tau=0.05", "CCI tau=0.1", "CCI tau=0.2", "GLC")

# ---- paths -------------------------------------------------------------------
area_path <- here("src", "area_0p25_validdomain_km2.nc")

outdir <- here("analysis", "results", "paper_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

unm_path <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_yearamp_0p25.nc", var)
)

path_cci <- function(tau) {
  here("output", tau, "eval", sprintf("trend_%s_CCI", var), sprintf("%s_yearamp_0p25.nc", var))
}
path_glc <- function() {
  here("output", tau_glc, "eval", sprintf("trend_%s_GLC", var), sprintf("%s_yearamp_0p25.nc", var))
}
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

# ---- load + compute ----------------------------------------------------------
rows <- list()

r_unm <- rast(unm_path)
z_unm <- zonal_wmean_latbands(time_mean(r_unm), area, band_deg = band_deg) |>
  as_tibble() |>
  rename(mean_yearamp = value) |>
  mutate(scenario = "Unmasked")
rows[[length(rows) + 1]] <- z_unm

for (tau in taus_cci) {
  r_cci <- rast(path_cci(tau))
  rows[[length(rows) + 1]] <-
    zonal_wmean_latbands(time_mean(r_cci), area, band_deg = band_deg) |>
    as_tibble() |>
    rename(mean_yearamp = value) |>
    mutate(scenario = sprintf("CCI %s", gsub("tau_", "tau=", tau)))
}

r_glc <- rast(path_glc())
rows[[length(rows) + 1]] <-
  zonal_wmean_latbands(time_mean(r_glc), area, band_deg = band_deg) |>
  as_tibble() |>
  rename(mean_yearamp = value) |>
  mutate(scenario = "GLC")

zonal_tbl <- bind_rows(rows) |>
  mutate(scenario = factor(scenario, levels = scenario_levels))
# ---- write table -------------------------------------------------------------
out_csv <- file.path(
  outdir,
  "zonal_yearamp_timeMean_LAI_all_masks_tau_0.1.csv"
)
write_csv(mutate(zonal_tbl, across(where(is.double), ~ round(.x, 4))), out_csv)

# ---- plot --------------------------------------------------------------------
z_abs <- zonal_tbl |>
  filter(is.finite(mean_yearamp))

col_abs <- c(
  "Unmasked" = "grey20",
  "CCI tau=0.05" = "#1b9e77",
  "CCI tau=0.1" = "#d95f02",
  "CCI tau=0.2" = "#7570b3",
  "GLC" = "#386cb0"
)

ylab_abs <- if (toupper(var) == "LAI") {
  "Mean seasonal amplitude (LAI)"
} else {
  "Mean seasonal amplitude (fAPAR)"
}
p_abs <- ggplot(z_abs, aes(lat_band, mean_yearamp, colour = scenario)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_line(linewidth = 0.85, na.rm = TRUE) +
  scale_x_continuous(
    limits = c(-90, 90),
    breaks = seq(-90, 90, by = 10),
    labels = lat_labels
  ) +
  scale_colour_manual(values = col_abs) +
  labs(x = NULL, y = ylab_abs, colour = NULL) +
  theme_pub(include_legend = TRUE, include_strip = TRUE)

p <- p_abs +
  labs(
    title = sprintf("%s: zonal seasonal amplitude across masking scenarios", var),
    subtitle = sprintf("Area-weighted zonal means (band = %d deg).", band_deg),
    x = "Latitude"
  )

out_png <- file.path(
  outdir,
  "zonal_yearamp_timeMean_LAI_all_masks_tau_0.1.png"
)
out_pdf <- file.path(
  outdir,
  "zonal_yearamp_timeMean_LAI_all_masks_tau_0.1.pdf"
)

ggsave(out_png, p, width = 10, height = 4.8, dpi = 320)
ggsave(out_pdf, p, width = 10, height = 4.8)
