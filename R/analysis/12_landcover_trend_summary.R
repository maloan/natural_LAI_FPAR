# ==============================================================================
# 12_landcover_trend_summary.R
# From yearly 0.25° majority LC maps -> stable LC (mode across years),
# then area-weighted trend summaries by LC class.
#
# Methodology:
#   1. Load annual 0.25° land-cover majority maps (1992-2020)
#   2. Merge LC subclasses to aggregated classes (IPCC level 1/2)
#   3. Compute stable class as mode of time series (with optional stability threshold)
#   4. Zonal statistics: area, mean trend per class (unmasked vs masked)
#   5. Bootstrap CIs for each class (n_boot=1000, 95% conf)
#   6. Output: paper-grade table + bar plot with retained area labels
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(here)
  library(ggplot2)
  library(tidyr)
  library(readr)
})

source(here("R", "helpers", "cli_args.R"))
source(here("R", "helpers", "bootstrap_ci.R"))

terraOptions(progress = 0, memfrac = 0.6, todisk = TRUE)

# ---- config ------------------------------------------------------------------
default_cfg <- list(
  tau = "tau_0.1",
  mask = "CCI",
  var = "LAI",
  metric = "yearmean",
  use_relative = FALSE,
  lc_year_start = 1992L,
  lc_year_end = 2020L,
  stability_min = NA_real_ # set numeric value (e.g., 0.6) to enable
)

cfg <- parse_cli_args(default_cfg)

tau <- as.character(cfg$tau)
mask <- as.character(cfg$mask)
var <- as.character(cfg$var)
metric <- as.character(cfg$metric)
use_relative <- isTRUE(cfg$use_relative)

if (is.na(cfg$lc_year_start) || is.na(cfg$lc_year_end) || cfg$lc_year_start > cfg$lc_year_end) {
  stop("Invalid lc year bounds: lc_year_start must be <= lc_year_end and both numeric")
}
lc_years <- cfg$lc_year_start:cfg$lc_year_end

stability_min <- if (is.na(cfg$stability_min)) NULL else cfg$stability_min

lc_yearly_dir <- here("analysis", "tmp", "lc025_majority_yearly")
out_dir <- here("analysis", "results", "lc_trends", tau)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
paper_fig_dir <- here("analysis", "results", "figures", "summaries")
dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Auto-generate landcover majority maps if missing
# ==============================================================================

lc_files_expected <- file.path(lc_yearly_dir, sprintf("lc025_majority_%d.tif", lc_years))
lc_files_missing <- lc_files_expected[!file.exists(lc_files_expected)]

if (length(lc_files_missing) > 0) {
  cat(sprintf(
    "\nNOTE: %d/%d landcover majority maps missing. Generating from ESACCI data...\n\n",
    length(lc_files_missing), length(lc_files_expected)
  ))

  lc_maker_script <- here("R", "12_make_lc025_majority.R")

  if (!file.exists(lc_maker_script)) {
    stop(
      sprintf(
        "Landcover preprocessing script not found: %s\nPlease run R/12_make_lc025_majority.R manually.",
        lc_maker_script
      ),
      call. = FALSE
    )
  }

  # Run the preprocessing script
  result <- system(sprintf("Rscript %s", lc_maker_script))

  if (result != 0) {
    stop(
      "Landcover preprocessing failed. Please check R/12_make_lc025_majority.R output.",
      call. = FALSE
    )
  }

  cat("\nLandcover preprocessing complete. Continuing with analysis...\n\n")
}

ref025 <- rast(here("src", "ref_0p25.nc"))
area <- rast(here("src", "area_0p25_validdomain_km2.nc"))[[1]]

# ---- helper ------------------------------------------------------------------
merge_to_parent <- c(
  # cropland + mosaics -> 10
  `11` = 10, `12` = 10, `20` = 10, `30` = 10, `40` = 10,

  # broadleaf group: 61/62 -> 60, then 60 -> 50
  `61` = 50, `62` = 50,
  `60` = 50,

  # needleleaf group: 71/72 -> 70; 81/82 -> 80; then 80 -> 70
  `71` = 70, `72` = 70,
  `81` = 70, `82` = 70,
  `80` = 70,

  # mosaic herbaceous/tree-shrub -> 100
  `110` = 100,

  # shrubland subclasses -> 120
  `121` = 120, `122` = 120,

  # sparse subclasses -> 150
  `151` = 150, `152` = 150, `153` = 150,

  # flooded -> 180
  `160` = 180, `170` = 180,

  # bare subclasses -> 200
  `201` = 200, `202` = 200
)

stable_mode <- function(v, min_frac = NULL) {
  v <- v[!is.na(v)]
  n <- length(v)
  if (!n) {
    return(NA_integer_)
  }
  ux <- unique(v)
  cnt <- tabulate(match(v, ux))
  k <- which.max(cnt)
  frac <- cnt[k] / n
  if (!is.null(min_frac) && frac < min_frac) {
    return(NA_integer_)
  }
  as.integer(ux[k])
}

align_to_ref <- function(r, ref, method = "near") {
  if (compareGeom(r, ref, stopOnError = FALSE)) {
    return(r)
  }
  resample(r, ref, method = method)
}

trend_files <- function(use_relative) {
  suf <- if (use_relative) "trend_relative_peryear" else "trend_slope_peryear"
  list(
    unm = here(
      "analysis", "unmasked", "0p25",
      sprintf("%s_georef_%s_%s_0p25.nc", var, metric, suf)
    ),
    msk = here(
      "output", tau, "eval", sprintf("trend_%s_%s", var, mask),
      sprintf("%s_%s_%s_0p25.nc", var, metric, suf)
    )
  )
}

# ---- (1) stable LC map -------------------------------------------------------
lc_files <- file.path(lc_yearly_dir, sprintf("lc025_majority_%d.tif", lc_years))

missing_files <- lc_files[!file.exists(lc_files)]
if (length(missing_files) > 0) {
  stop(
    sprintf(
      "Missing %d landcover files in %s.\nExpected files:\n  %s\n\nFirst missing: %s",
      length(missing_files),
      lc_yearly_dir,
      paste(basename(lc_files), collapse = "\n  "),
      basename(missing_files[1])
    ),
    call. = FALSE
  )
}

lc_stack <- rast(lc_files)
lc_stack <- align_to_ref(lc_stack, ref025, method = "near") |> mask(area)
# Merge subclasses on the 0.25° yearly maps
rcl <- cbind(
  from = as.integer(names(merge_to_parent)),
  to = as.integer(merge_to_parent)
)
lc_stack <- classify(lc_stack, rcl = rcl)
merge_tag <- "T1T2_v2"
lc_cache <- file.path(
  out_dir,
  sprintf(
    "lc025_stable_%s_%d-%d_min%s.tif",
    merge_tag, min(lc_years), max(lc_years),
    ifelse(is.null(stability_min), "none", gsub("\\.", "p", stability_min))
  )
)

lc_stable <- app(lc_stack, \(v) stable_mode(v, min_frac = stability_min),
  filename = lc_cache, overwrite = TRUE,
  wopt = list(datatype = "INT2S", gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6"))
)
names(lc_stable) <- "lc_id"


# ---- (2) read trends ---------------------------------------------------------
tf <- trend_files(use_relative)
stopifnot(file.exists(tf$unm), file.exists(tf$msk))

r_unm <- rast(tf$unm)[[1]] |> align_to_ref(ref025, method = "bilinear")
r_msk <- rast(tf$msk)[[1]] |> align_to_ref(ref025, method = "bilinear")

scale_factor <- if (use_relative) 100 else 1
suffix <- if (use_relative) "rel" else "abs"
unit_label <- if (use_relative) "% yr-1" else sprintf("%s yr-1", var)

# ---- (3) ONE zonal() call (fast) + bootstrap CIs ----
w_unm <- area
w_msk <- ifel(is.na(r_msk), NA_real_, area)
w_out <- ifel(!is.na(r_unm) & is.na(r_msk), area, NA_real_)

num_unm <- r_unm * w_unm
num_msk <- r_msk * w_msk
num_out <- r_unm * w_out

x <- c(num_unm, num_msk, num_out, w_unm, w_msk, w_out)
names(x) <- c("num_unm", "num_msk", "num_out", "den_unm_km2", "den_msk_km2", "den_out_km2")

zone <- zonal(x, lc_stable, fun = "sum", na.rm = TRUE) |> as.data.frame()
names(zone)[1] <- "lc_id"

# Compute bootstrap CIs for unmasked and masked trends
r_unm_vals <- terra::values(r_unm[[1]], dataframe = FALSE)
r_msk_vals <- terra::values(r_msk[[1]], dataframe = FALSE)
w_vals <- terra::values(w_unm[[1]], dataframe = FALSE)
z_vals <- terra::values(lc_stable[[1]], dataframe = FALSE)

ci_unm <- bootstrap_ci_by_class(r_unm_vals * scale_factor, z_vals, w_vals, n_boot = 1000L) |>
  rename(
    n_unm = n_pixels,
    mean_unmasked_boot = mean_est,
    ci_unm_lower = ci_lower,
    ci_unm_upper = ci_upper,
    ci_unm_width = ci_width
  )

ci_msk <- bootstrap_ci_by_class(r_msk_vals * scale_factor, z_vals, w_vals, n_boot = 1000L) |>
  rename(
    n_msk = n_pixels,
    mean_masked_boot = mean_est,
    ci_msk_lower = ci_lower,
    ci_msk_upper = ci_upper,
    ci_msk_width = ci_width
  )

lc_tab <- zone |>
  mutate(
    mean_unmasked   = scale_factor * (num_unm / den_unm_km2),
    mean_masked     = scale_factor * (num_msk / den_msk_km2),
    mean_masked_out = scale_factor * (num_out / den_out_km2),
    area_unm_mkm2   = den_unm_km2 / 1e6,
    area_msk_mkm2   = den_msk_km2 / 1e6,
    area_out_mkm2   = den_out_km2 / 1e6,
    frac_retained   = ifelse(den_unm_km2 > 0, den_msk_km2 / den_unm_km2, NA_real_)
  ) |>
  left_join(ci_unm, by = c("lc_id" = "class")) |>
  left_join(ci_msk, by = c("lc_id" = "class")) |>
  mutate(across(starts_with("mean_"), ~ ifelse(is.finite(.x), .x, NA_real_))) |>
  arrange(desc(area_unm_mkm2))

out_csv_full <- file.path(out_dir, sprintf("lc_class_%s_%s_%s_%s.csv", var, metric, mask, suffix))
write_csv(lc_tab, out_csv_full)

# ---- (4) plot (no legend work; just numeric classes) -------------------------
lc_legend_merged <- tibble::tribble(
  ~lc_id, ~lc_name,
  10L, "Cropland (incl. mosaics)",
  50L, "Tree cover, broadleaved",
  70L, "Tree cover, needleleaved",
  90L, "Tree cover, mixed",
  100L, "Mosaic tree/shrub and herbaceous",
  120L, "Shrubland",
  130L, "Grassland",
  140L, "Lichens and mosses",
  150L, "Sparse vegetation (<15%)",
  180L, "Flooded vegetation / wetlands",
  190L, "Urban areas",
  200L, "Bare areas",
  210L, "Water bodies",
  220L, "Permanent snow and ice"
)
lc_tab
plot_tab <- lc_tab |>
  filter(!is.na(lc_id)) |>
  left_join(lc_legend_merged, by = "lc_id") |>
  mutate(lc_name = coalesce(lc_name, paste0("Class ", lc_id))) |>
  arrange(desc(mean_unmasked))

plot_df <- plot_tab |>
  select(
    lc_name, mean_unmasked, mean_masked, ci_unm_lower, ci_unm_upper,
    ci_msk_lower, ci_msk_upper, frac_retained
  ) |>
  pivot_longer(
    cols = c(mean_unmasked, mean_masked),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    ci_lower = if_else(metric == "mean_unmasked", ci_unm_lower, ci_msk_lower),
    ci_upper = if_else(metric == "mean_unmasked", ci_unm_upper, ci_msk_upper),
    metric = factor(metric,
      levels = c("mean_unmasked", "mean_masked"),
      labels = c("Unmasked", "Masked")
    ),
    lc_name = factor(lc_name, levels = rev(unique(plot_tab$lc_name)))
  )

pd <- position_dodge(width = 0.7)
finite_vals <- plot_df$value[is.finite(plot_df$value)]
if (!length(finite_vals)) {
  yoff <- 0
} else {
  yspan <- diff(range(finite_vals))
  yoff <- 0.03 * ifelse(is.finite(yspan) && yspan > 0, yspan, max(abs(finite_vals)))
}

area_unm_tot <- sum(lc_tab$area_unm_mkm2, na.rm = TRUE)
area_msk_tot <- sum(lc_tab$area_msk_mkm2, na.rm = TRUE)

paper_tab <- plot_tab |>
  transmute(
    lc_id,
    lc_name,
    area_unmasked_mkm2 = area_unm_mkm2,
    area_masked_mkm2 = area_msk_mkm2,
    retained_pct = 100 * frac_retained,
    trend_unmasked = mean_unmasked,
    trend_unmasked_ci_lower = ci_unm_lower,
    trend_unmasked_ci_upper = ci_unm_upper,
    trend_masked = mean_masked,
    trend_masked_ci_lower = ci_msk_lower,
    trend_masked_ci_upper = ci_msk_upper,
    trend_delta = mean_masked - mean_unmasked,
    sig_unmasked = sig_flag.x,
    sig_masked = sig_flag.y
  ) |>
  mutate(across(where(is.double), ~ round(.x, 4)))

out_csv_paper <- file.path(
  paper_fig_dir,
  sprintf("table_lc_trend_summary_%s_%s_%s_%s_main.csv", var, metric, mask, tau)
)
write_csv(paper_tab, out_csv_paper)

base_labs <- labs(
  y = unit_label, fill = NULL,
  title = sprintf("Land-cover class trends (%s, %s)", var, metric),
  subtitle = sprintf(
    "%s mask at %s; stable classes from yearly 0.25 deg majority (%d-%d). Labels = retained area (%%)",
    mask, tau, min(lc_years), max(lc_years)
  ),
  caption = sprintf(
    "Total area: unmasked %.1f Mkm2, masked %.1f Mkm2",
    area_unm_tot, area_msk_tot
  )
)

if (!nrow(plot_df) || !length(finite_vals)) {
  p_bar <- ggplot() +
    annotate("text", x = 1, y = 1, label = "No finite class trends available") +
    xlim(0.5, 1.5) +
    ylim(0.5, 1.5) +
    theme_bw(base_size = 10) +
    base_labs +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold")
    )
} else {
  p_bar <- ggplot(plot_df, aes(lc_name, value, fill = metric)) +
    # Error bars with CI
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.15, position = pd, linewidth = 0.4, na.rm = TRUE,
      alpha = 0.7
    ) +
    geom_col(position = pd, width = 0.65, na.rm = TRUE) +
    geom_text(
      data = filter(plot_df, metric == "Masked"),
      aes(
        y = value + sign(value) * yoff,
        label = sprintf("%.0f%%", 100 * frac_retained)
      ),
      position = pd, size = 2.6, na.rm = TRUE
    ) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    scale_fill_manual(values = c("Unmasked" = "grey40", "Masked" = "black")) +
    base_labs +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "top",
      plot.margin = margin(5.5, 25, 5.5, 5.5),
      plot.title = element_text(face = "bold")
    )
}

out_png <- file.path(out_dir, sprintf("plot_lc_bar_%s_%s_%s_%s.png", var, metric, mask, suffix))
out_pdf <- file.path(out_dir, sprintf("plot_lc_bar_%s_%s_%s_%s.pdf", var, metric, mask, suffix))

paper_png <- file.path(paper_fig_dir, sprintf("lc_class_trend_summary_%s_%s_%s_%s_main.png", var, metric, mask, tau))
paper_pdf <- file.path(paper_fig_dir, sprintf("lc_class_trend_summary_%s_%s_%s_%s_main.pdf", var, metric, mask, tau))

ggsave(out_png, p_bar, width = 8.4, height = 6.2, dpi = 320)
ggsave(out_pdf, p_bar, width = 8.4, height = 6.2)
ggsave(paper_png, p_bar, width = 8.4, height = 6.2, dpi = 320)
ggsave(paper_pdf, p_bar, width = 8.4, height = 6.2)

cat("Wrote:\n")
