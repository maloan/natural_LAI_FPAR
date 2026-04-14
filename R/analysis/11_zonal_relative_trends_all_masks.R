# ==============================================================================
# 11_zonal_relative_trends_all_masks.R
# Figure 4 components: zonal relative trends for annual mean and annual max
# All scenarios overlaid (unmasked, CCI tau sweep, GLC)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(here)
})

# ---- config ------------------------------------------------------------------
var <- "LAI"
metrics <- c("yearmean", "yearmax")
taus_cci <- c("tau_0.05", "tau_0.1", "tau_0.2")
tau_glc <- "tau_0.1"

outdir <- here("analysis", "results", "paper_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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

check_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing ", label, ": ", path)
  }
}

check_file_exists(area_path, "area raster")
area <- rast(area_path)[[1]]

path_unmasked <- function(metric) {
  here(
    "analysis", "unmasked", "0p25",
    sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
}
path_cci <- function(metric, tau) {
  here(
    "output", tau, "eval", sprintf("trend_%s_CCI", var),
    sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
}
path_glc <- function(metric) {
  here(
    "output", tau_glc, "eval", sprintf("trend_%s_GLC", var),
    sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
}

# ---- zonal mean helper -------------------------------------------------------
lat_labels <- function(x) {
  ifelse(
    x == 0, "0°",
    ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N"))
  )
}
lat_breaks_30 <- seq(-90, 90, by = 30)

# Returns lat band midpoints, zonal mean (%/yr), contributing area (km2)
zonal_wmean_lat <- function(r, area) {
  if (nlyr(r) > 1) {
    r <- r[[1]]
  }

  compareGeom(r, area, stopOnError = TRUE)

  lat <- init(area, "y")
  zone <- floor((lat + 90) / 1) + 1

  ok <- is.finite(r) & is.finite(area) & (area > 0)
  num <- ifel(ok, r * area, NA)
  den <- ifel(ok, area, NA)

  z_num <- zonal(num, zone, "sum", na.rm = TRUE)
  names(z_num) <- c("zone", "num_sum")
  z_den <- zonal(den, zone, "sum", na.rm = TRUE)
  names(z_den) <- c("zone", "den_sum")

  z <- merge(z_num, z_den, by = "zone", all = TRUE)
  z$lat <- -90 + (z$zone - 0.5) * 1

  z$reltrend_pct_per_year <- 100 * (z$num_sum / z$den_sum)

  as_tibble(z) |>
    transmute(
      lat_band = .data$lat,
      reltrend_pct_per_year = ifelse(is.finite(.data$reltrend_pct_per_year),
        .data$reltrend_pct_per_year, NA_real_
      ),
      area_km2 = ifelse(is.finite(.data$den_sum), .data$den_sum, NA_real_)
    ) |>
    arrange(.data$lat_band)
}
# ---- compute -----------------------------------------------------------------
rows <- list()

for (met in metrics) {
  f_unmasked <- path_unmasked(met)
  check_file_exists(f_unmasked, sprintf("unmasked relative-trend raster (%s)", met))
  r_unmasked <- rast(f_unmasked)[[1]]
  rows[[length(rows) + 1]] <- zonal_wmean_lat(r_unmasked, area) |>
    mutate(metric = met, scenario = "Unmasked")

  for (tau in taus_cci) {
    f_cci <- path_cci(met, tau)
    check_file_exists(f_cci, sprintf("CCI relative-trend raster (%s, %s)", met, tau))
    r_cci <- rast(f_cci)[[1]]
    rows[[length(rows) + 1]] <- zonal_wmean_lat(r_cci, area) |>
      mutate(metric = met, scenario = sprintf("CCI %s", gsub("tau_", "tau=", tau)))
  }

  f_glc <- path_glc(met)
  check_file_exists(f_glc, sprintf("GLC relative-trend raster (%s, %s)", met, tau_glc))
  r_glc <- rast(f_glc)[[1]]
  rows[[length(rows) + 1]] <- zonal_wmean_lat(r_glc, area) |>
    mutate(metric = met, scenario = "GLC")
}

df <- bind_rows(rows) |>
  mutate(
    metric = factor(metric, levels = c("yearmean", "yearmax"), labels = c("Annual mean", "Annual maximum")),
    scenario = factor(scenario, levels = c("Unmasked", "CCI tau=0.05", "CCI tau=0.1", "CCI tau=0.2", "GLC"))
  )

# ---- write table -------------------------------------------------------------
write_csv(df, file.path(outdir, sprintf("zonal_relative_trends_all_masks_%s.csv", tau_glc)))

# ---- plot --------------------------------------------------------------------
col_map <- c(
  "Unmasked" = "grey20",
  "CCI tau=0.05" = "#1b9e77",
  "CCI tau=0.1" = "#d95f02",
  "CCI tau=0.2" = "#7570b3",
  "GLC" = "#386cb0"
)

base_plot <- function(df_plot, ttl) {
  ggplot(df_plot, aes(.data$lat_band, .data$reltrend_pct_per_year, colour = .data$scenario)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_line(linewidth = 0.65, na.rm = TRUE) +
    scale_colour_manual(values = col_map, drop = FALSE) +
    scale_x_continuous(
      limits = c(-90, 90),
      breaks = lat_breaks_30,
      labels = lat_labels
    ) +
    labs(
      x = "Latitude",
      y = expression(paste("% yr"^{
        -1
      })),
      colour = NULL,
      title = ttl
    ) +
    theme_pub()
}

p_mean <- base_plot(
  df |> filter(metric == "Annual mean"),
  "Zonal relative trends: annual mean"
)

p_max <- base_plot(
  df |> filter(metric == "Annual maximum"),
  "Zonal relative trends: annual maximum"
)

ggsave(file.path(outdir, "zonal_relative_trends_yearmean_all_masks.png"), p_mean, width = 10, height = 4.7, dpi = 320)
ggsave(file.path(outdir, "zonal_relative_trends_yearmean_all_masks.pdf"), p_mean, width = 10, height = 4.7)
ggsave(file.path(outdir, "zonal_relative_trends_yearmax_all_masks.png"), p_max, width = 10, height = 4.7, dpi = 320)
ggsave(file.path(outdir, "zonal_relative_trends_yearmax_all_masks.pdf"), p_max, width = 10, height = 4.7)
