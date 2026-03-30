# ==============================================================================
# global_timeseries_masks.R
# Global annual time series (unmasked + CCI tau sweep + GLC) — 0.25° only
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(here)
  library(scales)
  library(scico)
})

# ---- config ------------------------------------------------------------------
root <- Sys.getenv("root", here::here())
outdir <- Sys.getenv("outdir", file.path(root, "analysis", "results", "global_timeseries_masks"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

vars <- strsplit(Sys.getenv("vars", "LAI"), ",")[[1]] |> trimws()
metrics <- strsplit(Sys.getenv("metrics", "yearmean,yearmax"), ",")[[1]] |> trimws()
taus_cci <- strsplit(Sys.getenv("taus_cci", "tau_0.05,tau_0.1,tau_0.2"), ",")[[1]] |> trimws()
tau_glc <- Sys.getenv("tau_glc", "tau_0.1")

in_geo <- file.path(root, "analysis", "unmasked", "0p25")

# ---- weights -----------------------------------------------------------------
area025_path <- file.path(root, "src", "area_0p25_validdomain_km2.nc")
area005_path <- file.path(root, "src", "area_0p05_validdomain_km2.nc")
stopifnot(file.exists(area025_path), file.exists(area005_path))

area_abi_025 <- rast(area025_path)[[1]] # post-abiotic land base at 0.25°
area_abi_005 <- rast(area005_path)[[1]] # post-abiotic land base at 0.05°

abiotic_land_km2 <- global(ifel(is.finite(area_abi_025) & (area_abi_025 > 0),
                                area_abi_025, NA), "sum", na.rm = TRUE)[1, 1] |>
  as.numeric()

# ---- helpers -----------------------------------------------------------------
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

metric_labels <- c(yearmean = "Annual mean", yearmax = "Annual maximum")
stopifnot(all(metrics %in% names(metric_labels)))

make_label <- function(domain_type, mask, tau) {
  if (domain_type == "unmasked") {
    return("Unmasked (post-abiotic)")
  }
  if (mask == "GLC") {
    return("Natural-only (GLC-based)")
  }
  sprintf("Natural-only (CCI-based) (tau=%s)", sub("^tau_", "", tau))
}

label_levels <- function(taus_cci) {
  c(
    "Unmasked (post-abiotic)",
    sprintf("Natural-only (CCI-based) (tau=%s)", sub("^tau_", "", taus_cci)),
    "Natural-only (GLC-based)"
  )
}

# ---- mask paths --------------------------------------------------------------
mask_path_cci <- function(tau) {
  tau_num <- as.numeric(sub("^tau_", "", tau))
  tau_token <- sprintf("tau0p%02d", round(tau_num * 100))
  file.path(
    root, "output", tau, "masks", "mask_cci",
    sprintf("mask_used_frac_fused_%s_k3_1992-2020_0p05.tif", tau_token)
  )
}

mask_path_glc <- function(tau_glc) {
  file.path(root, "output", tau_glc, "masks", "mask_glc", "mask_used_ge3_1992-2020_0p05.tif")
}

# Build remaining-area weights at 0.25° from a 0.05° (possibly fractional) mask.
# mask005 is interpreted as "excluded fraction" (0..1), so remaining = (1 - mask).
eff_area_025_from_mask005 <- function(mask005_path, area_abi_005, area_abi_025_template) {
  stopifnot(file.exists(mask005_path))
  m005 <- rast(mask005_path)[[1]]
  compareGeom(m005, area_abi_005, stopOnError = TRUE)

  ok005 <- is.finite(area_abi_005) & (area_abi_005 > 0) & is.finite(m005)
  area_eff005 <- ifel(ok005, (1 - m005) * area_abi_005, NA)

  area_eff025 <- aggregate(area_eff005, fact = 5, fun = "sum", na.rm = TRUE)
  compareGeom(area_eff025, area_abi_025_template, stopOnError = TRUE)

  area_eff025
}

# Global mean time series using provided weights area_weight (km²).
# Returns:
#   value   : area-weighted mean of x over valid x (within area_weight)
#   area_km2: effective support area (sum of area_weight over valid x)
global_mean_series <- function(r, area_weight, year0 = 1982L) {
  yrs <- year0:(year0 + nlyr(r) - 1L)
  out_val <- numeric(nlyr(r))
  out_den <- numeric(nlyr(r))

  for (i in seq_len(nlyr(r))) {
    x <- r[[i]]
    ok <- is.finite(x) & is.finite(area_weight) & (area_weight > 0)

    den <- global(ifel(ok, area_weight, NA), "sum", na.rm = TRUE)[1, 1] |> as.numeric()
    if (!is.finite(den) || den <= 0) {
      out_val[i] <- NA_real_
      out_den[i] <- NA_real_
      next
    }

    num <- global(ifel(ok, x * area_weight, NA), "sum", na.rm = TRUE)[1, 1] |> as.numeric()
    out_val[i] <- num / den
    out_den[i] <- den
  }

  tibble(year = yrs, value = out_val, area_km2 = out_den)
}

# ---- precompute effective remaining-area weights (static, per mask) ----------
area_eff_list <- list()

# CCI τ sweep
for (tau in taus_cci) {
  fmask <- mask_path_cci(tau)
  stopifnot(file.exists(fmask))
  area_eff_list[[paste0("CCI__", tau)]] <- eff_area_025_from_mask005(fmask, area_abi_005, area_abi_025)
}

# GLC (single)
fmask <- mask_path_glc(tau_glc)
stopifnot(file.exists(fmask))
area_eff_list[[paste0("GLC__", tau_glc)]] <- eff_area_025_from_mask005(fmask, area_abi_005, area_abi_025)


static_remaining_km2 <- function(area_eff025) {
  global(ifel(is.finite(area_eff025) & (area_eff025 > 0), area_eff025, NA), "sum", na.rm = TRUE)[1, 1] |>
    as.numeric()
}

# ---- main --------------------------------------------------------------------
for (var in vars) {
  rows <- list()

  # --- unmasked (weights = full post-abiotic area at 0.25°) -------------------
  for (met in metrics) {
    f <- file.path(in_geo, sprintf("%s_georef_%s_0p25.nc", var, met))
    stopifnot(file.exists(f))
    rows[[length(rows) + 1]] <-
      global_mean_series(rast(f), area_abi_025) |>
      mutate(
        domain_type = "unmasked",
        mask = NA_character_,
        tau = NA_character_,
        var = var,
        metric = met,
        source_file = f,
        area_static_km2 = abiotic_land_km2,
        area_static_remaining_km2 = abiotic_land_km2, # unmasked: remaining = abiotic base
        area_static_excluded_km2 = 0 # unmasked: excluded = 0
      )
  }

  # --- CCI masked (weights = remaining-area weights derived from mask) --------
  for (tau in taus_cci) {
    area_eff025 <- area_eff_list[[paste0("CCI__", tau)]]
    area_rem_static <- static_remaining_km2(area_eff025)
    area_excl_static <- abiotic_land_km2 - area_rem_static

    for (met in metrics) {
      f <- file.path(root, "output", tau, "eval", sprintf("trend_%s_CCI", var), sprintf("%s_%s_0p25.nc", var, met))
      stopifnot(file.exists(f))

      rows[[length(rows) + 1]] <-
        global_mean_series(rast(f), area_eff025) |>
        mutate(
          domain_type = "masked",
          mask = "CCI",
          tau = tau,
          var = var,
          metric = met,
          source_file = f,
          area_static_km2 = abiotic_land_km2, # always report the common abiotic base
          area_static_remaining_km2 = area_rem_static, # static remaining after mask (km²)
          area_static_excluded_km2 = area_excl_static # static excluded by mask (km²)
        )
    }
  }

  # --- GLC masked --------------------------------------------------------------
  area_eff025 <- area_eff_list[[paste0("GLC__", tau_glc)]]
  area_rem_static <- static_remaining_km2(area_eff025)
  area_excl_static <- abiotic_land_km2 - area_rem_static

  for (met in metrics) {
    f <- file.path(root, "output", tau_glc, "eval", sprintf("trend_%s_GLC", var), sprintf("%s_%s_0p25.nc", var, met))
    stopifnot(file.exists(f))

    rows[[length(rows) + 1]] <-
      global_mean_series(rast(f), area_eff025) |>
      mutate(
        domain_type = "masked",
        mask = "GLC",
        tau = tau_glc,
        var = var,
        metric = met,
        source_file = f,
        area_static_km2 = abiotic_land_km2,
        area_static_remaining_km2 = area_rem_static,
        area_static_excluded_km2 = area_excl_static
      )
  }

  df <- bind_rows(rows) |>
    mutate(
      label = mapply(make_label, domain_type, mask, tau),
      label = factor(label, levels = label_levels(taus_cci)),
      mask_group = case_when(
        domain_type == "unmasked" ~ "unmasked",
        mask == "GLC" ~ "GLC",
        mask == "CCI" ~ "CCI",
        TRUE ~ NA_character_
      ),
      mask_group = factor(mask_group, levels = c("unmasked", "CCI", "GLC")),
      metric = factor(metric, levels = names(metric_labels), labels = unname(metric_labels)),

      # Effective support shortfall relative to the *static remaining footprint*
      # (this is the quantity you should interpret as missing observations)
      area_missing_km2 = pmax(0, area_static_remaining_km2 - area_km2),

      # Convenience: static masked pct (theoretical footprint)
      footprint_masked_pct = 100 * area_static_excluded_km2 / area_static_km2
    ) |>
    arrange(label, metric, year)

  # ---- colours / linetypes ---------------------------------------------------
  labs <- levels(df$label)
  is_unm <- labs == "Unmasked (post-abiotic)"
  is_glc <- labs == "Natural-only (GLC-based)"
  is_cci <- grepl("^Natural-only \\(CCI-based\\) \\(tau=", labs)

  cols <- setNames(rep(NA_character_, length(labs)), labs)
  cols[is_unm] <- "grey15"
  cols[is_glc] <- "#7570b3"
  if (sum(is_cci) > 0) cols[is_cci] <- scico(sum(is_cci), palette = "batlow")

  linetypes <- c(unmasked = "solid", CCI = "solid", GLC = "dashed")

  sizes <- setNames(rep(0.60, length(labs)), labs)
  sizes["Unmasked (post-abiotic)"] <- 0.90
  sizes["Natural-only (GLC-based)"] <- 0.75

  # ---- plot ------------------------------------------------------------------
  p <- ggplot(df, aes(year, value, colour = label, group = label)) +
    geom_line(aes(linewidth = label, linetype = mask_group), na.rm = TRUE) +
    scale_linewidth_manual(values = sizes, guide = "none") +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_linetype_manual(values = linetypes, drop = FALSE) +
    facet_wrap(~metric, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(1982, max(df$year, na.rm = TRUE), by = 5)) +
    labs(
      x = "Year",
      y = sprintf("Absolute %s", var),
      title = sprintf("%s: global annual time series (masked vs unmasked)", var),
      subtitle = paste0(
        "Unmasked series uses the post-abiotic land base as weights.\n",
        "Masked series use mask-derived remaining-area weights (static footprint), and the support diagnostic is\n",
        "(static remaining area − effective observed support area)."
      ),
      colour = NULL,
      linetype = "Series type"
    ) +
    theme_pub() +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))

  out_base <- file.path(outdir, sprintf("%s_global_timeseries_yearmean_yearmax_masks", var))

  ggsave(paste0(out_base, ".png"), p, width = 7.6, height = 6.6, dpi = 400)
  write_csv(df, paste0(out_base, ".csv"))

  # Effective support per year (km²)
  write_csv(
    df |>
      select(year, metric, label, mask_group, area_km2) |>
      arrange(metric, label, year),
    paste0(out_base, "_area_km2.csv")
  )

  # Effective vs static diagnostics (km²)
  write_csv(
    df |>
      select(
        year, metric, label, mask_group,
        area_km2,
        area_static_km2,
        area_static_remaining_km2,
        area_static_excluded_km2,
        footprint_masked_pct,
        area_missing_km2
      ) |>
      arrange(metric, label, year),
    paste0(out_base, "_area_effective_vs_static_km2.csv")
  )

  # Summary by series
  write_csv(
    df |>
      group_by(metric, label) |>
      summarise(
        area_static_km2 = first(area_static_km2),
        area_static_remaining_km2 = first(area_static_remaining_km2),
        area_static_excluded_km2 = first(area_static_excluded_km2),
        footprint_masked_pct = first(footprint_masked_pct),
        area_km2_min = min(area_km2, na.rm = TRUE),
        area_km2_mean = mean(area_km2, na.rm = TRUE),
        area_km2_max = max(area_km2, na.rm = TRUE),
        area_missing_km2_mean = mean(area_missing_km2, na.rm = TRUE),
        .groups = "drop"
      ),
    paste0(out_base, "_area_effective_vs_static_summary.csv")
  )
}
