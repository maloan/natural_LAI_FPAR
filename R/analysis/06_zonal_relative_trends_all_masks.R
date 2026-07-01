# ==============================================================================
# 06_zonal_relative_trends_all_masks.R — zonal relative trends for annual mean
# and annual max
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(here)
})

source(here("R", "helpers", "weighted_means.R"))
source(here("R", "helpers", "plotting.R"))
source(here("R", "helpers", "bootstrap_ci.R"))
source(here("R", "helpers", "io.R"))

# config
var <- "LAI"
metrics <- c("yearmean", "yearmax")
taus_cci <- c("tau_0.05", "tau_0.1", "tau_0.2")
tau_glc <- "tau_0.1"
outdir_fig <- here("analysis", "results", "figures", "summaries")
outdir_tbl <- here("analysis", "results", "tables", "zonal")
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tbl, recursive = TRUE, showWarnings = FALSE)

# paths
area_path <- here("src", "area_0p25_validdomain_km2.nc")
check_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing ", label, ": ", path)
  }
}
check_file_exists(area_path, "area raster")
area <- rast(area_path)[[1]]
block_id <- make_block_id(area, block_size_deg = 5)

zonal_wmean_latbands_ci <- function(r,
                                    area,
                                    block_id,
                                    band_deg = 1L,
                                    scale_factor = 1,
                                    n_boot = 500L,
                                    conf = 0.95) {
  lat <- terra::init(area, "y")
  lat_vals <- terra::values(lat, dataframe = FALSE)
  r_vals <- terra::values(r, dataframe = FALSE)
  area_vals <- terra::values(area, dataframe = FALSE)
  lat_min <- seq(-90, 90 - band_deg, by = band_deg)
  lat_max <- lat_min + band_deg
  results <- vector("list", length(lat_min))
  for (i in seq_along(lat_min)) {
    in_band <- lat_vals >= lat_min[i] & lat_vals < lat_max[i]
    r_band <- r_vals[in_band] * scale_factor
    area_band <- area_vals[in_band]
    block_band <- block_id[in_band]
    ok <- is.finite(r_band) &
      is.finite(area_band) & area_band > 0 & !is.na(block_band)
    if (sum(ok) < 2) {
      results[[i]] <- tibble::tibble(
        lat_band = lat_min[i] + band_deg / 2,
        value = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        n_pixels = sum(ok)
      )
      next
    }
    ci <- bootstrap_ci_global(
      x = r_band[ok],
      w = area_band[ok],
      block_id = block_band[ok],
      n_boot = n_boot,
      conf = conf
    )
    results[[i]] <- tibble::tibble(
      lat_band = lat_min[i] + band_deg / 2,
      value = ci$mean,
      ci_lower = ci$lower,
      ci_upper = ci$upper,
      n_pixels = ci$n_eff
    )
  }
  dplyr::bind_rows(results)
}
rows <- list()
for (met in metrics) {
  f_unmasked <- analysis_raster_path(var, met, "unmasked", kind = "trend_relative")
  check_file_exists(
    f_unmasked,
    sprintf("unmasked relative-trend raster (%s)", met)
  )
  r_unmasked <- load_checked_raster(f_unmasked, area, label = "Unmasked", first_layer = TRUE)
  rows[[length(rows) + 1]] <- zonal_wmean_latbands_ci(
    r_unmasked,
    area,
    block_id,
    band_deg = 1L,
    scale_factor = 100,
    n_boot = 500L
  ) |>
    as_tibble() |>
    rename(reltrend_pct_per_year = value) |>
    arrange(.data$lat_band) |>
    mutate(metric = met, scenario = "Unmasked")

  for (tau in taus_cci) {
    f_cci <- analysis_raster_path(var, met, "CCI", run_tag = tau, kind = "trend_relative")
    check_file_exists(
      f_cci,
      sprintf("CCI relative-trend raster (%s, %s)", met, tau)
    )
    r_cci <- load_checked_raster(f_cci,
      area,
      label = paste("CCI", tau),
      first_layer = TRUE
    )
    rows[[length(rows) + 1]] <- zonal_wmean_latbands_ci(
      r_cci,
      area,
      block_id,
      band_deg = 1L,
      scale_factor = 100,
      n_boot = 500L
    ) |>
      as_tibble() |>
      rename(reltrend_pct_per_year = value) |>
      arrange(.data$lat_band) |>
      mutate(metric = met, scenario = sprintf("CCI %s", gsub("tau_", "tau=", tau)))
  }

  f_glc <- analysis_raster_path(var, met, "GLC", run_tag = tau_glc, kind = "trend_relative")
  check_file_exists(
    f_glc,
    sprintf("GLC relative-trend raster (%s, %s)", met, tau_glc)
  )
  r_glc <- load_checked_raster(f_glc, area, label = "GLC", first_layer = TRUE)
  rows[[length(rows) + 1]] <- zonal_wmean_latbands_ci(
    r_glc,
    area,
    block_id,
    band_deg = 1L,
    scale_factor = 100,
    n_boot = 500L
  ) |>
    as_tibble() |>
    rename(reltrend_pct_per_year = value) |>
    arrange(.data$lat_band) |>
    mutate(metric = met, scenario = "GLC")
}
df <- bind_rows(rows) |>
  mutate(
    metric = factor(
      metric,
      levels = c("yearmean", "yearmax"),
      labels = c("Annual mean", "Annual maximum")
    ),
    scenario = factor(
      scenario,
      levels = c("Unmasked", "CCI tau=0.05", "CCI tau=0.1", "CCI tau=0.2", "GLC")
    )
  )

# write table
write_csv(round_numeric(df, 5), file.path(
  outdir_tbl,
  sprintf("zonal_relative_trends_all_masks_%s.csv", tau_glc)
))

p_mean <- plot_zonal_reltrend(
  df |> filter(metric == "Annual mean"),
  "Zonal relative trends: annual mean"
)
p_max <- plot_zonal_reltrend(
  df |> filter(metric == "Annual maximum"),
  "Zonal relative trends: annual maximum"
)
ggsave(
  file.path(outdir_fig, "zonal_relative_trends_yearmean_all_masks.png"),
  p_mean,
  width = 10,
  height = 4.7,
  dpi = 320
)
ggsave(
  file.path(outdir_fig, "zonal_relative_trends_yearmean_all_masks.pdf"),
  p_mean,
  width = 10,
  height = 4.7
)
ggsave(
  file.path(outdir_fig, "zonal_relative_trends_yearmax_all_masks.png"),
  p_max,
  width = 10,
  height = 4.7,
  dpi = 320
)
ggsave(
  file.path(outdir_fig, "zonal_relative_trends_yearmax_all_masks.pdf"),
  p_max,
  width = 10,
  height = 4.7
)
