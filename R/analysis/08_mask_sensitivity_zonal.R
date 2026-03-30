# ==============================================================================
# 08_mask_sensitivity_tau.R
# Sensitivity of headline greening to mask strictness (tau)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(here)
})

# ---- config ------------------------------------------------------------------
metric <- Sys.getenv("metric", "yearmean")
taus <- strsplit(("tau_0.05,tau_0.1,tau_0.2"), ",")[[1]]
vars <- "LAI"

outdir <- here("analysis", "results", "mask_sensitivity_tau")
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
      legend.position  = "right",
      legend.box       = "vertical",
      legend.text      = element_text(size = 9)
    )
}
# ---- paths -------------------------------------------------------------------
area025 <- here("src", "area_0p25_validdomain_km2.nc")
area <- rast(area025)[[1]]
valid_dom <- is.finite(area) & (area > 0)

f_unmasked <- function(var) {
  here(
    "analysis", "unmasked", "0p25",
    sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
}
f_masked_cci <- function(tau, var) {
  here(
    "output", tau, "eval", sprintf("trend_%s_CCI", var),
    sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
}
# ---- helpers -----------------------------------------------------------------
lat_labels <- function(x) {
  ifelse(
    x == 0, "0°",
    ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N"))
  )
}
lat_breaks_30 <- seq(-90, 90, by = 30)

wmean_global <- function(x, w) {
  ok <- valid_dom & is.finite(x) & is.finite(w) & (w > 0)
  den <- global(ifel(ok, w, NA), "sum", na.rm = TRUE)[1, 1]
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  num <- global(ifel(ok, x * w, NA), "sum", na.rm = TRUE)[1, 1]
  as.numeric(num / den)
}

zonal_wmean <- function(x, w) {
  # assumes x and w aligned
  lat <- init(w, "y")
  zone <- floor((lat + 90) / 1) + 1

  ok <- valid_dom & is.finite(x) & is.finite(w) & (w > 0)
  num <- ifel(ok, x * w, NA)
  den <- ifel(ok, w, NA)

  s_num <- zonal(num, zone, "sum", na.rm = TRUE)
  names(s_num) <- c("zone", "num_sum")
  s_den <- zonal(den, zone, "sum", na.rm = TRUE)
  names(s_den) <- c("zone", "den_sum")

  out <- merge(s_num, s_den, by = "zone", all = TRUE)

  tibble(
    lat = -90 + (out$zone - 0.5),
    reltrend_pct_per_year = 100 * (out$num_sum / out$den_sum),
    area_km2 = out$den_sum
  ) |>
    mutate(
      reltrend_pct_per_year = ifelse(is.finite(.data$reltrend_pct_per_year),
        .data$reltrend_pct_per_year, NA_real_
      ),
      area_km2 = ifelse(is.finite(.data$area_km2), .data$area_km2, NA_real_)
    ) |>
    arrange(lat)
}

# ---- palette -----------------------------------------------------------------
# Unmasked in grey; τ series in sequential blues (light -> dark).
pal_series <- c(
  "unmasked" = "grey20",
  "tau_0.05" = "#a6cee3",
  "tau_0.1"  = "#1f78b4",
  "tau_0.2"  = "#08306b"
)

lab_series <- c("unmasked" = "Unmasked (post-abiotic domain)")
lab_series <- c(lab_series, setNames(gsub("^tau_", "CCI τ = ", taus), taus))

# ---- compute -----------------------------------------------------------------
zonal_rows <- list()
global_rows <- list()

for (var in vars) {
  fu <- f_unmasked(var)
  rU <- rast(fu)[[1]]
  # Global (unmasked)
  global_rows[[length(global_rows) + 1]] <- tibble(
    variable = var, metric = metric, series = "unmasked",
    tau = NA_character_, tau_value = NA_real_,
    global_reltrend_pct_per_year = 100 * wmean_global(rU, area),
    area_km2 = as.numeric(
      global(ifel(valid_dom, area, NA), "sum", na.rm = TRUE)[1, 1]
    )
  )

  # Zonal (unmasked)
  zonal_rows[[length(zonal_rows) + 1]] <- zonal_wmean(rU, area) |>
    mutate(
      variable = var, metric = metric, series = "unmasked",
      tau = NA_character_, tau_value = NA_real_
    )

  for (tau in taus) {
    fc <- f_masked_cci(tau, var)
    rC <- rast(fc)[[1]]

    global_rows[[length(global_rows) + 1]] <- tibble(
      variable = var, metric = metric, series = tau,
      tau = tau, tau_value = as.numeric(sub("^tau_", "", tau)),
      global_reltrend_pct_per_year = 100 * wmean_global(rC, area),
      area_km2 = as.numeric(
        global(
          ifel(valid_dom & is.finite(rC), area, NA), "sum",
          na.rm = TRUE
        )[1, 1]
      )
    )

    zonal_rows[[length(zonal_rows) + 1]] <-
      zonal_wmean(rC, area) |>
      mutate(
        variable = var, metric = metric, series = tau,
        tau = tau, tau_value = as.numeric(sub("^tau_", "", tau))
      )
  }
}

zon <- bind_rows(zonal_rows) |>
  mutate(
    variable = factor(variable, levels = c("LAI", "FPAR")),
    series = factor(series, levels = c("unmasked", taus))
  )

glob <- bind_rows(global_rows) |>
  mutate(
    variable = factor(variable, levels = c("LAI", "FPAR")),
    series = factor(series, levels = c("unmasked", taus))
  )

# ---- write CSVs --------------------------------------------------------------
write_csv(zon, file.path(
  outdir, sprintf("mask_sensitivity_zonal_%s.csv", metric)
))
write_csv(glob, file.path(
  outdir, sprintf("mask_sensitivity_global_%s.csv", metric)
))

# Paper-facing global table (LAI default; keep tidy even if vars has more)
tab_global <- glob |>
  filter(metric == metric) |>
  transmute(
    Variable = as.character(variable),
    Series = as.character(series),
    `Global mean relative trend (% yr-1)` = round(
      global_reltrend_pct_per_year, 3
    ),
    `Effective area (km2)` = round(area_km2, 0)
  ) |>
  arrange(Variable, factor(Series, levels = c("unmasked", taus)))

write_csv(tab_global, file.path(
  outdir, sprintf("table_mask_sensitivity_global_%s.csv", metric)
))

# ---- plot: global sensitivity ------------------------------------------------
pg <- ggplot(
  glob,
  aes(x = series, y = global_reltrend_pct_per_year, colour = series)
) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_point(size = 2.4, na.rm = TRUE) +
  geom_line(aes(group = variable), linewidth = 0.6, na.rm = TRUE) +
  facet_grid(. ~ variable, scales = "free_y") +
  scale_colour_manual(values = pal_series, labels = lab_series, drop = FALSE) +
  labs(
    x = NULL,
    y = "Global mean relative trend (% yr-1)",
    colour = "Series",
    title = sprintf(
      "Global sensitivity to CCI mask strictness (metric: %s)", metric
    )
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(
  file.path(outdir, sprintf("mask_sensitivity_global_%s.png", metric)),
  pg,
  width = 8.8, height = 3.6, dpi = 300
)

# ---- plot: zonal sensitivity -------------------------------------------------
pz <- ggplot(
  zon,
  aes(lat, reltrend_pct_per_year, colour = series, group = series)
) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  facet_grid(. ~ variable, scales = "free_y") +
  scale_colour_manual(values = pal_series, labels = lab_series, drop = FALSE) +
  scale_x_continuous(
    limits = c(-90, 90),
    breaks = lat_breaks_30,
    labels = lat_labels
  ) +
  labs(
    x = "Latitude",
    y = "Zonal mean relative trend (% yr-1)",
    colour = "Series",
    title = sprintf(
      "Latitudinal sensitivity to CCI mask strictness (metric: %s; lat bins: %g°)",
      metric, 1
    ),
    subtitle = "Note: effective area per latitude band is recorded in the CSV (area_km2)."
  ) +
  theme_pub()


ggsave(
  file.path(outdir, sprintf("mask_sensitivity_zonal_%s.png", metric)),
  pz,
  width = 10.5, height = 4.2, dpi = 300
)
