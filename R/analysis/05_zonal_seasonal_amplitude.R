# ==============================================================================
# 05_zonal_seasonal_amplitude.R — Zonal mean seasonal amplitude (absolute
# values)
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
source(here("R", "helpers", "io.R"))
# config
var <- "LAI"
alphas_cci <- c("alpha_0.05", "alpha_0.1", "alpha_0.2")
alpha_glc <- "alpha_0.1"
band_deg <- 1L
scenario_levels <- c("Unmasked", "CCI alpha=0.05", "CCI alpha=0.1", "CCI alpha=0.2", "GLC")
# paths
area_path <- here("src", "area_0p25_validdomain_km2.nc")
outdir_fig <- here("analysis", "results", "figures", "summaries")
outdir_tbl <- here("analysis", "results", "tables", "zonal")
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tbl, recursive = TRUE, showWarnings = FALSE)

area <- rast(area_path)[[1]]

# time-mean area raster that may have multiple layers (annual time series)
time_mean <- function(r) {
  if (nlyr(r) == 1) {
    return(r[[1]])
  }
  terra::mean(r, na.rm = TRUE)
}
# load data and compute zonal summaries
rows <- list()
r_unm <- load_checked_raster(analysis_raster_path(var, "yearamp", "unmasked", kind = "metric"),
                             area,
                             label = "Unmasked")
z_unm <- zonal_wmean_latbands(time_mean(r_unm), area, band_deg = band_deg) |>
  as_tibble() |>
  rename(mean_yearamp = value) |>
  mutate(scenario = "Unmasked")
rows[[length(rows) + 1]] <- z_unm
for (alpha in alphas_cci) {
  r_cci <- load_checked_raster(
    analysis_raster_path(var, "yearamp", "CCI", run_tag = alpha, kind = "metric"),
    area,
    label = alpha
  )
  rows[[length(rows) + 1]] <-
    zonal_wmean_latbands(time_mean(r_cci), area, band_deg = band_deg) |>
    as_tibble() |>
    rename(mean_yearamp = value) |>
    mutate(scenario = sprintf("CCI %s", gsub("alpha_", "alpha=", alpha)))
}
r_glc <- load_checked_raster(
  analysis_raster_path(var, "yearamp", "GLC", run_tag = alpha_glc, kind = "metric"),
  area,
  label = "GLC"
)
rows[[length(rows) + 1]] <-
  zonal_wmean_latbands(time_mean(r_glc), area, band_deg = band_deg) |>
  as_tibble() |>
  rename(mean_yearamp = value) |>
  mutate(scenario = "GLC")
zonal_tbl <- bind_rows(rows) |>
  mutate(scenario = factor(scenario, levels = scenario_levels))
# plot
z_abs <- zonal_tbl |>
  filter(is.finite(mean_yearamp))
p <- plot_seasonal_amplitude(z_abs)
out_png <- file.path(outdir_fig,
                     sprintf("zonal_yearamp_timeMean_%s_all_masks_alpha_0.1.png", var))
out_pdf <- file.path(outdir_fig,
                     sprintf("zonal_yearamp_timeMean_%s_all_masks_alpha_0.1.pdf", var))
# write output
out_csv <- file.path(outdir_tbl,
                     sprintf("zonal_yearamp_timeMean_%s_all_masks_alpha_0.1.csv", var))
write_csv(round_numeric(zonal_tbl, 5), out_csv)
ggsave(out_png,
       p,
       width = 10,
       height = 4.8,
       dpi = 320)
ggsave(out_pdf, p, width = 10, height = 4.8)
