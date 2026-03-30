# ==============================================================================
# 08_zonal_relative_trends.R
# Zonal means of relative greening (independent per product + Δ)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(here)
  library(patchwork)
})

# ---- config ------------------------------------------------------------------
tau <- "tau_0.1"
vars <- c("LAI") # or c("LAI","FPAR")
metrics <- c("yearmean", "yearmax")
masktag <- "CCI" # or "GLC"

outdir <- here("analysis", "results", "zonal_relative_trends", tau)
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
area <- rast(area_path)[[1]]

path_unmasked <- function(var, metric) {
  here(
    "analysis", "unmasked", "0p25",
    sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
}
path_masked <- function(var, metric, masktag) {
  here(
    "output", tau, "eval",
    sprintf("trend_%s_%s", var, masktag),
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
      lat_band = lat,
      reltrend_pct_per_year = ifelse(is.finite(.data$reltrend_pct_per_year),
        .data$reltrend_pct_per_year, NA_real_
      ),
      area_km2 = ifelse(is.finite(.data$den_sum), .data$den_sum, NA_real_)
    ) |>
    arrange(.data$lat_band)
}
# ---- compute -----------------------------------------------------------------
rows <- list()

for (var in vars) {
  for (met in metrics) {
    fU <- path_unmasked(var, met)
    fM <- path_masked(var, met, masktag)

    rU <- rast(fU)[[1]]
    rM <- rast(fM)[[1]]

    zU <- zonal_wmean_lat(rU, area) |>
      mutate(case = "unmasked")
    zM <- zonal_wmean_lat(rM, area) |>
      mutate(case = paste0("masked_", masktag))

    z <- bind_rows(zU, zM) |>
      mutate(variable = var, metric = met)

    zD <- zU |>
      select(lat_band, reltrend_pct_per_year, area_km2) |>
      rename(rel_unm = reltrend_pct_per_year, area_unm_km2 = area_km2) |>
      left_join(
        zM |>
          select(lat_band, reltrend_pct_per_year, area_km2) |>
          rename(rel_msk = reltrend_pct_per_year, area_msk_km2 = area_km2),
        by = "lat_band"
      ) |>
      transmute(
        lat_band,
        reltrend_pct_per_year = rel_msk - rel_unm,
        area_unm_km2,
        area_msk_km2,
        case = "delta"
      ) |>
      mutate(variable = var, metric = met)

    # For consistency, carry area_km2 as NA for delta row
    z <- z |>
      mutate(area_unm_km2 = NA_real_, area_msk_km2 = NA_real_) |>
      bind_rows(zD |>
                  mutate(area_km2 = NA_real_))
  }
  rows[[length(rows) + 1]] <- z
}

df <- bind_rows(rows) |>
  mutate(
    variable = factor(variable, levels = c("LAI", "FPAR")),
    metric   = factor(metric, levels = metrics),
    case     = factor(case, levels = c("unmasked", paste0("masked_", masktag), "delta"))
  )

# ---- write table -------------------------------------------------------------
write_csv(df, file.path(outdir, sprintf("zonal_relative_trends_%s_%s.csv", masktag, tau)))

# ---- plot --------------------------------------------------------------------
lab_main <- c(
  "Unmasked (post-abiotic)",
  sprintf("Natural-only (%s-based)", masktag)
)
names(lab_main) <- c("unmasked", paste0("masked_", masktag))
col_main <- c("grey20", "#7570b3")
names(col_main) <- c("unmasked", paste0("masked_", masktag))

df_main <- df |>
  filter(case %in% names(lab_main)) |>
  mutate(case = droplevels(case))

p_abs <- ggplot(df_main, aes(lat_band, reltrend_pct_per_year, colour = case)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_line(linewidth = 0.65, na.rm = TRUE) +
  facet_grid(metric ~ variable, scales = "fixed") +
  scale_colour_manual(values = col_main, labels = lab_main, drop = FALSE) +
  scale_x_continuous(
    limits = c(-90, 90),
    breaks = lat_breaks_30,
    labels = lat_labels
  ) +
  labs(x = NULL, y = expression(paste("% yr"^{
    -1
  })), colour = NULL) +
  theme_pub()

df_d <- df |> filter(case == "delta")
p_del <- ggplot(df_d, aes(lat_band, reltrend_pct_per_year)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
  geom_area(fill = "#D55E00", alpha = 0.25, na.rm = TRUE) +
  geom_line(colour = "#D55E00", linewidth = 0.65, na.rm = TRUE) +
  facet_grid(metric ~ variable, scales = "fixed") +
  scale_x_continuous(
    limits = c(-90, 90),
    breaks = lat_breaks_30,
    labels = lat_labels
  ) +
  labs(x = "Latitude", y = expression(paste(Delta, " (% yr"^{
    -1
  }, ")"))) +
  theme_pub() +
  theme(legend.position = "none")


p <- p_abs / p_del + plot_layout(heights = c(1.15, 1.0)) +
  plot_annotation(
    title = sprintf(
                    "%s: zonal relative greening before/after masking (%s; %s)",
                    paste(vars, collapse = ","), masktag, tau),
    subtitle =
      "Zonal means computed independently per product; Δ is the difference of zonal means (band = 1°)."
  )


ggsave(
  file.path(outdir, sprintf("zonal_relative_trends_%s_%s.png", masktag, tau)),
  p,
  width = 10, height = if (ADD_DELTA) 8.0 else 6.0, dpi = 300
)
