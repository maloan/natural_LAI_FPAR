# ==============================================================================
# 07_dropped_region_signal.R
# Greening signal removed by land-use masking
#
# Definition used here:
#   dropped = cells that are valid in the unmasked trend (post-abiotic domain)
#             but are NA in the masked trend output.
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(rnaturalearth)
  library(sf)
  library(scales)
  library(here)
  library(scico)
})

# ---- config ------------------------------------------------------------------
tau <- Sys.getenv("tau", "tau_0.1")
vars <- "LAI"
metric <- "yearmean"
masks <- "CCI"

outdir <- here("analysis", "results", "dropped_region", tau)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- style -------------------------------------------------------------------
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

path_unmasked <- function(var) {
  here(
    "analysis", "unmasked", "0p25",
    sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
}

path_masked <- function(var, masktag) {
  here(
    "output", tau, "eval",
    sprintf("trend_%s_%s", var, masktag),
    sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
}

# ---- helpers -----------------------------------------------------------------
lat_labels <- function(x) {
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N")))
}
lon_labels <- function(x) {
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°W"), paste0(x, "°E")))
}

lat_breaks_30 <- seq(-90, 90, by = 30)
lon_breaks_60 <- seq(-180, 180, by = 30)

wmean <- function(x, m) {
  ok <- m & valid_dom & is.finite(x) & is.finite(area) & (area > 0)
  den <- global(ifel(ok, area, NA), "sum", na.rm = TRUE)[1, 1]
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  num <- global(ifel(ok, x * area, NA), "sum", na.rm = TRUE)[1, 1]
  as.numeric(num / den)
}

robust_zlim <- function(r, q = 0.98, fallback = 5, min_lim = 1, min_n = 200L) {
  v <- values(r, mat = FALSE)
  v <- v[is.finite(v)]
  if (length(v) < min_n) {
    return(fallback)
  }
  z <- as.numeric(stats::quantile(abs(v), probs = q, na.rm = TRUE))
  if (!is.finite(z)) fallback else max(min_lim, z)
}

coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

plot_map <- function(r_pct, title, out_png, zlim) {
  df <- as.data.frame(r_pct, xy = TRUE, na.rm = FALSE)
  names(df) <- c("lon", "lat", "val")

  p <- ggplot(df) +
    geom_raster(aes(x = .data$lon, y = .data$lat, fill = .data$val)) +
    geom_sf(data = coast, linewidth = 0.2, colour = "black", fill = NA, inherit.aes = FALSE) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    scale_x_continuous(breaks = lon_breaks_60, labels = lon_labels) +
    scale_y_continuous(breaks = lat_breaks_30, labels = lat_labels) +
    scale_fill_gradientn(
      colours = scico::scico(256, palette = "bam", direction = 1),
      limits = c(-zlim, zlim),
      oob = squish,
      na.value = "transparent",
      name = "% yr-1"
    ) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_pub()

  ggsave(out_png, p, width = 10.5, height = 5.2, dpi = 300)
}

zonal_wmean <- function(x, m) {
  lat <- init(area, "y")
  zone <- floor((lat + 90) / 1) + 1

  ok <- m & valid_dom & is.finite(x) & is.finite(area) & (area > 0)
  num <- ifel(ok, x * area, NA)
  den <- ifel(ok, area, NA)

  s_num <- zonal(num, zone, "sum", na.rm = TRUE)
  names(s_num) <- c("zone", "num_sum")
  s_den <- zonal(den, zone, "sum", na.rm = TRUE)
  names(s_den) <- c("zone", "den_sum")

  out <- merge(s_num, s_den, by = "zone", all = TRUE) |>
    tibble::as_tibble() |>
    mutate(
      num_sum = dplyr::coalesce(.data$num_sum, 0),
      den_sum = dplyr::coalesce(.data$den_sum, 0)
    )

  rel <- ifelse(out$den_sum > 0, 100 * (out$num_sum / out$den_sum), NA_real_)

  tibble(
    lat = -90 + (out$zone - 0.5),
    reltrend_pct_per_year = rel,
    area_km2 = out$den_sum
  ) |>
    arrange(lat)
}


plot_zonal <- function(df, title, out_png) {
  p <- ggplot(df, aes(.data$lat, .data$reltrend_pct_per_year)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
    geom_line(linewidth = 0.8, colour = "black", na.rm = TRUE) +
    scale_x_continuous(breaks = lat_breaks_30, labels = lat_labels) +
    labs(title = title, x = "Latitude (deg)", y = "Dropped-region rel. trend (% yr-1)") +
    theme_pub()

  ggsave(out_png, p, width = 7.5, height = 4.0, dpi = 300)
}

# ---- compute -----------------------------------------------------------------
rows <- list()

area_land_km2 <- as.numeric(global(ifel(valid_dom, area, NA), "sum", na.rm = TRUE)[1, 1])

for (var in vars) {
  fU <- path_unmasked(var)
  rU <- rast(fU)[[1]]

  area_unm_valid_km2 <- as.numeric(
    global(ifel(is.finite(rU) & valid_dom, area, NA), "sum", na.rm = TRUE)[1, 1]
  )

  # Stable zlim reference from unmasked field over post-abiotic domain
  ru_pct_dom <- ifel(valid_dom & is.finite(rU), 100 * rU, NA_real_)
  zlim_ref <- robust_zlim(ru_pct_dom, q = 0.98, fallback = 5)

  for (masktag in masks) {
    fM <- path_masked(var, masktag)
    rM <- rast(fM)[[1]]

    # --- dropped region (THIS is the definition you asked for) ----------------
    dropped <- valid_dom & is.finite(rU) & is.na(rM)

    # areas
    area_drop_km2 <- as.numeric(global(ifel(dropped, area, NA), "sum", na.rm = TRUE)[1, 1])

    drop_pct_land <- if (is.finite(area_land_km2) && area_land_km2 > 0) {
      100 * area_drop_km2 / area_land_km2
    } else {
      NA_real_
    }

    drop_pct_unmvalid <- if (is.finite(area_unm_valid_km2) && area_unm_valid_km2 > 0) {
      100 * area_drop_km2 / area_unm_valid_km2
    } else {
      NA_real_
    }

    # mean signal in dropped region (based on unmasked field)
    mu_drop_pctyr <- 100 * wmean(rU, dropped)

    message(sprintf(
      "[%s %s %s] dropped_area=%.0f km2 (%.2f%% of post-abiotic land; %.2f%% of unmasked-valid); mean=%.3f %%/yr",
      var, metric, masktag,
      area_drop_km2, drop_pct_land, drop_pct_unmvalid, mu_drop_pctyr
    ))

    rows[[length(rows) + 1]] <- tibble(
      tau = tau, variable = var, metric = metric, mask = masktag,
      denom_postabiotic_land_km2 = area_land_km2,
      denom_unmasked_valid_km2 = area_unm_valid_km2,
      dropped_area_km2 = area_drop_km2,
      dropped_area_pct_of_land = drop_pct_land,
      dropped_area_pct_of_unmaskedvalid = drop_pct_unmvalid,
      dropped_mean_reltrend_pct_per_year = mu_drop_pctyr
    )

    # Map: unmasked reltrend restricted to dropped region (percent per year)
    r_drop_pct <- ifel(dropped, 100 * rU, NA_real_)
    map_png <- file.path(
      outdir,
      sprintf("map_dropped_reltrend_%s_%s_%s_%s.png", var, metric, masktag, tau)
    )
    plot_map(
      r_drop_pct,
      sprintf("%s %s — dropped region (%s): unmasked relative trend", var, metric, masktag),
      map_png,
      zlim_ref
    )

    # Zonal curve + per-band effective area
    z <- zonal_wmean(rU, dropped)
    z_csv <- file.path(
      outdir,
      sprintf("zonal_dropped_reltrend_%s_%s_%s_%s.csv", var, metric, masktag, tau)
    )
    z_png <- file.path(
      outdir,
      sprintf("zonal_dropped_reltrend_%s_%s_%s_%s.png", var, metric, masktag, tau)
    )

    write_csv(
      mutate(z, variable = var, metric = metric, mask = masktag, tau = tau),
      z_csv
    )
    z_plot <- dplyr::filter(z, area_km2 > 0)

    plot_zonal(
      z_plot,
      sprintf("%s %s — dropped region zonal mean (%s)", var, metric, masktag),
      z_png
    )
  }
}

summary_df <- bind_rows(rows)
write_csv(summary_df, file.path(outdir, sprintf("dropped_region_summary_%s.csv", tau)))

# ---- paper-facing table (yearmean only) --------------------------------------
tab_main <- summary_df |>
  filter(variable == "LAI", metric == "yearmean") |>
  transmute(
    tau = tau,
    Mask = mask,
    Variable = variable,
    `Dropped area (% of post-abiotic land)` = dropped_area_pct_of_land,
    `Mean relative trend in dropped region (% yr-1)` = dropped_mean_reltrend_pct_per_year
  ) |>
  mutate(across(where(is.double), ~ round(.x, 3))) |>
  arrange(Mask)

write_csv(
  tab_main,
  file.path(outdir, sprintf("table_dropped_region_LAI_yearmean_%s.csv", tau))
)
