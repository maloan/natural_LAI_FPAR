# ==============================================================================
# 06_zonal_relative_trends_all_masks.R
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

source(here("R", "helpers", "weighted_means.R"))
source(here("R", "helpers", "plotting.R"))

plot_theme_pub <- function(base_size = 12, include_legend = FALSE,
                           include_strip = FALSE) {
  theme_list <- list(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = base_size + 1, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = base_size - 1),
    axis.text = element_text(size = 9)
  )
  if (include_strip) {
    theme_list$strip.text <- element_text(size = 10, face = "bold")
  }
  if (include_legend) {
    theme_list$legend.position <- "bottom"
    theme_list$legend.box <- "vertical"
    theme_list$legend.text <- element_text(size = 9)
  }
  theme_bw(base_size = base_size) + theme(!!!theme_list)
}

theme_pub <- plot_theme_pub

# ---- Helper: Zonal weighted mean with bootstrap CIs ----
zonal_wmean_latbands_ci <- function(r, area, band_deg = 1L, scale_factor = 1,
                                    n_boot = 500L, conf = 0.95) {
  lat_min <- seq(-90, 90 - band_deg, by = band_deg)
  lat_max <- lat_min + band_deg

  results <- list()
  set.seed(42)

  for (i in seq_along(lat_min)) {
    # Extract latitude band pixels
    lat <- terra::init(area, "y")
    in_band <- terra::values(lat) >= lat_min[i] & terra::values(lat) < lat_max[i]

    r_band <- terra::values(r, dataframe = FALSE)[in_band]
    area_band <- terra::values(area, dataframe = FALSE)[in_band]

    ok <- is.finite(r_band) & is.finite(area_band) & area_band > 0
    if (sum(ok) == 0) {
      results[[i]] <- tibble(
        lat_band = lat_min[i] + band_deg / 2,
        value = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        n_pixels = 0
      )
      next
    }

    r_ok <- r_band[ok] * scale_factor
    w_ok <- area_band[ok]
    n <- sum(ok)

    # Point estimate
    wmean_est <- sum(r_ok * w_ok) / sum(w_ok)

    # Bootstrap resampling
    boot_means <- numeric(n_boot)
    for (j in seq_len(n_boot)) {
      boot_idx <- sample.int(n, size = n, replace = TRUE, prob = w_ok / sum(w_ok))
      boot_means[j] <- sum(r_ok[boot_idx] * w_ok[boot_idx]) / sum(w_ok)
    }

    # Percentile-based CI
    ci_alpha <- 1 - conf
    ci_lower <- quantile(boot_means, ci_alpha / 2, type = 7)
    ci_upper <- quantile(boot_means, 1 - ci_alpha / 2, type = 7)

    results[[i]] <- tibble(
      lat_band = lat_min[i] + band_deg / 2,
      value = wmean_est,
      ci_lower = as.numeric(ci_lower),
      ci_upper = as.numeric(ci_upper),
      n_pixels = n
    )
  }

  bind_rows(results)
}

# ---- config ------------------------------------------------------------------
var <- "LAI"
metrics <- c("yearmean", "yearmax")
taus_cci <- c("tau_0.05", "tau_0.1", "tau_0.2")
tau_glc <- "tau_0.1"

outdir <- here("analysis", "results", "figures", "summaries")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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

# ---- compute (with bootstrap CIs) ---------
rows <- list()

for (met in metrics) {
  f_unmasked <- path_unmasked(met)
  check_file_exists(f_unmasked, sprintf(
    "unmasked relative-trend raster (%s)",
    met
  ))
  r_unmasked <- rast(f_unmasked)[[1]]
  rows[[length(rows) + 1]] <- zonal_wmean_latbands_ci(r_unmasked,
    area,
    band_deg = 1L, scale_factor = 100, n_boot = 500L
  ) |>
    as_tibble() |>
    rename(reltrend_pct_per_year = value) |>
    arrange(.data$lat_band) |>
    mutate(metric = met, scenario = "Unmasked")

  for (tau in taus_cci) {
    f_cci <- path_cci(met, tau)
    check_file_exists(f_cci, sprintf(
      "CCI relative-trend raster (%s, %s)",
      met, tau
    ))
    r_cci <- rast(f_cci)[[1]]
    rows[[length(rows) + 1]] <- zonal_wmean_latbands_ci(r_cci,
      area,
      band_deg = 1L, scale_factor = 100, n_boot = 500L
    ) |>
      as_tibble() |>
      rename(reltrend_pct_per_year = value) |>
      arrange(.data$lat_band) |>
      mutate(metric = met, scenario = sprintf(
        "CCI %s",
        gsub("tau_", "tau=", tau)
      ))
  }

  f_glc <- path_glc(met)
  check_file_exists(f_glc, sprintf(
    "GLC relative-trend raster (%s, %s)",
    met, tau_glc
  ))
  r_glc <- rast(f_glc)[[1]]
  rows[[length(rows) + 1]] <- zonal_wmean_latbands_ci(r_glc, area,
    band_deg = 1L, scale_factor = 100, n_boot = 500L
  ) |>
    as_tibble() |>
    rename(reltrend_pct_per_year = value) |>
    arrange(.data$lat_band) |>
    mutate(metric = met, scenario = "GLC")
}

df <- bind_rows(rows) |>
  mutate(
    metric = factor(metric,
      levels = c("yearmean", "yearmax"),
      labels = c("Annual mean", "Annual maximum")
    ),
    scenario = factor(scenario,
      levels = c(
        "Unmasked", "CCI tau=0.05",
        "CCI tau=0.1", "CCI tau=0.2", "GLC"
      )
    )
  )

# ---- write table -------------------------------------------------------------
write_csv(df, file.path(
  outdir,
  sprintf("zonal_relative_trends_all_masks_%s.csv", tau_glc)
))

# ---- plot --------------------------------------------------------------------
col_map <- c(
  "Unmasked" = "grey20",
  "CCI tau=0.05" = "#1b9e77",
  "CCI tau=0.1" = "#d95f02",
  "CCI tau=0.2" = "#7570b3",
  "GLC" = "#386cb0"
)

base_plot <- function(df_plot, ttl) {
  ggplot(df_plot, aes(.data$lat_band, .data$reltrend_pct_per_year,
    colour = .data$scenario
  )) +
    # Confidence band (shaded ribbon)
    geom_ribbon(
      aes(ymin = .data$ci_lower, ymax = .data$ci_upper, fill = .data$scenario),
      alpha = 0.15, colour = NA, na.rm = TRUE
    ) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_line(linewidth = 0.65, na.rm = TRUE) +
    scale_colour_manual(values = col_map, drop = FALSE) +
    scale_fill_manual(values = col_map, drop = FALSE, guide = "none") +
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
      title = ttl,
      subtitle = "Lines: zonal mean; bands: 95% bootstrap CI"
    ) +
    plot_theme_pub(include_legend = TRUE, include_strip = TRUE)
}

p_mean <- base_plot(
  df |> filter(metric == "Annual mean"),
  "Zonal relative trends: annual mean"
)

p_max <- base_plot(
  df |> filter(metric == "Annual maximum"),
  "Zonal relative trends: annual maximum"
)

ggsave(file.path(outdir, "zonal_relative_trends_yearmean_all_masks.png"),
  p_mean,
  width = 10, height = 4.7, dpi = 320
)
ggsave(file.path(outdir, "zonal_relative_trends_yearmean_all_masks.pdf"),
  p_mean,
  width = 10, height = 4.7
)
ggsave(file.path(outdir, "zonal_relative_trends_yearmax_all_masks.png"),
  p_max,
  width = 10, height = 4.7, dpi = 320
)
ggsave(file.path(outdir, "zonal_relative_trends_yearmax_all_masks.pdf"),
  p_max,
  width = 10, height = 4.7
)
