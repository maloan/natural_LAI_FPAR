# ==============================================================================
# 12_a_landcover_trend_summary.R — Dominant land-cover class trend summaries
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(here)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
  library(readr)
  library(purrr)
})

source(here::here("R", "helpers", "cli_args.R"))
source(here::here("R", "helpers", "bootstrap_ci.R"))
source(here::here("R", "helpers", "plotting.R"))
source(here::here("R", "helpers", "io.R"))
terraOptions(
  progress = 0,
  memfrac = 0.6,
  todisk = TRUE
)

# config
default_cfg <- list(
  alpha = "alpha_0.1",
  mask = "CCI",
  var = "LAI",
  metric = "yearmean",
  use_relative = TRUE,
  lc_year_start = 1992L,
  lc_year_end = 2022L,
  n_boot = 1000L,
  conf = 0.95
)

cfg <- parse_cli_args(default_cfg)

alpha <- as.character(cfg$alpha)
mask <- as.character(cfg$mask)
var <- as.character(cfg$var)
metric <- as.character(cfg$metric)
use_relative <- isTRUE(cfg$use_relative)
n_boot <- as.integer(cfg$n_boot)
conf <- as.numeric(cfg$conf)

if (is.na(cfg$lc_year_start) ||
  is.na(cfg$lc_year_end) ||
  cfg$lc_year_start > cfg$lc_year_end) {
  stop("Invalid lc year bounds: lc_year_start must be <= lc_year_end and both numeric",
    call. = FALSE
  )
}

lc_years <- cfg$lc_year_start:cfg$lc_year_end

# paths
lc_class_dir <- here("analysis", "tmp", "lc025_dominant_class")
outdir_fig <- here("analysis", "results", "figures", "summaries")
outdir_tbl <- here("analysis", "results", "tables", "land_cover")

dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tbl, recursive = TRUE, showWarnings = FALSE)

ref025 <- terra::rast(here::here("src", "ref_0p25.nc"))
area_land <- rast(here::here("src", "area_0p25_validdomain_km2.nc"))[[1]]
area <- load_summary_area(template = area_land)
area_vals <- terra::values(area, dataframe = FALSE)
land_vals <- terra::values(area_land, dataframe = FALSE)
valid_domain_cells <- which(is.finite(land_vals) & land_vals > 0)
block_id <- make_block_id(area, block_size_deg = 5)

trend_files <- function(use_relative) {
  suf <- if (use_relative) {
    "trend_relative_percent_peryear"
  } else {
    "trend_slope_peryear"
  }

  list(
    unm = here(
      "analysis",
      "unmasked",
      "0p25",
      sprintf("%s_georef_%s_%s_0p25.nc", var, metric, suf)
    ),
    msk = here(
      "output",
      alpha,
      "eval",
      sprintf("trend_%s_%s", var, mask),
      sprintf("%s_%s_%s_0p25.nc", var, metric, suf)
    )
  )
}
# land-cover legend
lc_legend <- tibble::tribble(
  ~lc_id,
  ~lc_name,
  0L,
  "No Data",
  10L,
  "Cropland",
  50L,
  "Broadleaved Evergreen Forest",
  60L,
  "Broadleaved Deciduous Forest",
  70L,
  "Needleleaved Evergreen Forest",
  80L,
  "Needleleaved Deciduous Forest",
  90L,
  "Mixed Forest",
  100L,
  "Mosaic Tree/Shrub",
  120L,
  "Shrubland",
  130L,
  "Grassland",
  140L,
  "Lichen and Mosses",
  150L,
  "Sparse Vegetation",
  160L,
  "Flooded Tree Cover",
  180L,
  "Flooded Shrub/Herbaceous Cover",
  190L,
  "Urban Areas",
  200L,
  "Bare Areas",
  210L,
  "Water Bodies",
  220L,
  "Permanent Snow and Ice"
)

# Land-cover classes considered non-vegetated for plotting (exclude from
# figures)
nonveg_lc_ids <- c(
  200L,
  # Bare areas
  210L,
  # Water bodies
  220L # Permanent snow and ice
)

# Read the fixed dominant class. Land-cover fractions are used only for class
# assignment and are not statistical weights.
dominant_lc_file <- file.path(
  lc_class_dir,
  sprintf("lc025_dominant_mean_%d-%d.tif", min(lc_years), max(lc_years))
)
if (!file.exists(dominant_lc_file)) {
  stop(
    "Missing dominant land-cover classification. Run ",
    "R/analysis/12_0_landcover_dominant_class.R first:\n",
    dominant_lc_file
  )
}
dominant_lc <- rast(dominant_lc_file)[[1]]
stopifnot(compareGeom(dominant_lc, ref025, stopOnError = TRUE))
dominant_lc_values <- as.integer(terra::values(dominant_lc, dataframe = FALSE))
lc_ids <- sort(unique(dominant_lc_values[is.finite(dominant_lc_values)]))

# read trends
tf <- trend_files(use_relative)
r_unm <- rast(tf$unm)[[1]]
r_msk <- rast(tf$msk)[[1]]
r_unm <- align_to_template(r_unm, ref025, method = "bilinear")
r_msk <- align_to_template(r_msk, ref025, method = "bilinear")
r_unm_all <- terra::values(r_unm, dataframe = FALSE)
r_msk_all <- terra::values(r_msk, dataframe = FALSE)

if (use_relative) {
  scale_factor <- 1
  plot_scale_factor <- 1
  suffix <- "rel"
  unit_label <- "% yr-1"
} else {
  scale_factor <- 1
  plot_scale_factor <- 1000
  suffix <- "abs"
  unit_label <- expression("LAI trend (" ~ "×" ~ 10^-3 ~ m^2 ~ m^-2 ~ yr^-1 * ")")
}

# Each valid class member receives its fixed post-nonvegetated support area.
# Excluded-domain summaries remain restricted to cells without a masked trend.
summarise_dominant_class <- function(lc_id) {
  member <- dominant_lc_values == lc_id
  w_unm_cls <- ifelse(member & is.finite(r_unm_all), area_vals, NA_real_)
  w_msk_cls <- ifelse(member & is.finite(r_msk_all), area_vals, NA_real_)
  w_out_cls <- ifelse(member & is.finite(r_unm_all) & !is.finite(r_msk_all), area_vals, NA_real_)

  den_unm <- sum(w_unm_cls, na.rm = TRUE)
  den_msk <- sum(w_msk_cls, na.rm = TRUE)
  den_out_full_cells <- sum(w_out_cls, na.rm = TRUE)
  den_excluded <- max(den_unm - den_msk, 0)

  num_unm <- sum(r_unm_all * w_unm_cls, na.rm = TRUE)
  num_msk <- sum(r_msk_all * w_msk_cls, na.rm = TRUE)
  num_out_full_cells <- sum(r_unm_all * w_out_cls, na.rm = TRUE)

  tibble::tibble(
    lc_id = lc_id,
    mean_unmasked = scale_factor * safe_division(num_unm, den_unm),
    mean_masked = scale_factor * safe_division(num_msk, den_msk),
    mean_masked_out = scale_factor * safe_division(
      num_out_full_cells,
      den_out_full_cells
    ),
    area_unm_mkm2 = den_unm / 1e6,
    area_msk_mkm2 = den_msk / 1e6,
    area_out_mkm2 = den_excluded / 1e6,
    area_fully_excluded_cells_mkm2 = den_out_full_cells / 1e6,
    frac_retained = safe_division(den_msk, den_unm)
  )
}

lc_tab <- purrr::map_dfr(lc_ids, function(id) {
  summarise_dominant_class(id)
}) |>
  mutate(across(starts_with("mean_"), ~ ifelse(is.finite(.x), .x, NA_real_))) |>
  left_join(lc_legend, by = "lc_id") |>
  mutate(lc_name = coalesce(lc_name, paste0("Class ", lc_id))) |>
  arrange(desc(abs(mean_masked)))

# Subset shared bootstrap inputs once.
r_unm_vals <- r_unm_all[valid_domain_cells] *
  scale_factor
r_msk_vals <- r_msk_all[valid_domain_cells] *
  scale_factor
block_id_sub <- block_id[valid_domain_cells]

pb <- txtProgressBar(
  min = 0,
  max = length(lc_ids),
  style = 3
)

ci_list <- purrr::map2_dfr(seq_along(lc_ids), lc_ids, function(i, id) {
  setTxtProgressBar(pb, i)
  member_vals <- dominant_lc_values[valid_domain_cells] == id
  area_sub <- area_vals[valid_domain_cells]

  w_unm_cls <- ifelse(member_vals & is.finite(r_unm_vals), area_sub, NA_real_)
  w_msk_cls <- ifelse(member_vals & is.finite(r_msk_vals), area_sub, NA_real_)
  w_out_cls <- ifelse(
    member_vals & is.finite(r_unm_vals) & !is.finite(r_msk_vals),
    area_sub,
    NA_real_
  )

  ok_unm <- is.finite(r_unm_vals) &
    is.finite(w_unm_cls) & w_unm_cls > 0 & !is.na(block_id_sub)
  ci_unm <- bootstrap_ci_global(r_unm_vals[ok_unm], w_unm_cls[ok_unm], block_id_sub[ok_unm], n_boot, conf)
  ok_msk <- is.finite(r_msk_vals) &
    is.finite(w_msk_cls) & w_msk_cls > 0 & !is.na(block_id_sub)
  ci_msk <- bootstrap_ci_global(r_msk_vals[ok_msk], w_msk_cls[ok_msk], block_id_sub[ok_msk], n_boot, conf)

  ok_out <- is.finite(r_unm_vals) &
    is.finite(w_out_cls) & w_out_cls > 0 & !is.na(block_id_sub)
  ci_out <- bootstrap_ci_global(r_unm_vals[ok_out], w_out_cls[ok_out], block_id_sub[ok_out], n_boot, conf)

  difference_ci <- bootstrap_ci_difference(
    x_retained = r_msk_vals,
    w_retained = w_msk_cls,
    x_excluded = r_unm_vals,
    w_excluded = w_out_cls,
    block_id = block_id_sub,
    n_boot = n_boot,
    conf = conf
  )

  # Domain-only contrast: use the unmasked trend field in both groups and use
  # the mask only to classify 0.25-degree cells as retained or excluded.
  domain_difference_ci <- bootstrap_ci_difference(
    x_retained = r_unm_vals,
    w_retained = w_msk_cls,
    x_excluded = r_unm_vals,
    w_excluded = w_out_cls,
    block_id = block_id_sub,
    n_boot = n_boot,
    conf = conf
  )

  tibble::tibble(
    lc_id = id,
    ci_unm_lower = ci_unm$lower,
    ci_unm_upper = ci_unm$upper,
    ci_msk_lower = ci_msk$lower,
    ci_msk_upper = ci_msk$upper,
    ci_out_lower = ci_out$lower,
    ci_out_upper = ci_out$upper,
    n_unm = as.integer(ci_unm$n_eff),
    n_msk = as.integer(ci_msk$n_eff),
    retained_minus_excluded = difference_ci$estimate,
    difference_ci_lower = difference_ci$lower,
    difference_ci_upper = difference_ci$upper,
    difference_excludes_zero = difference_ci$excludes_zero,
    n_difference_blocks = as.integer(difference_ci$n_blocks),
    domain_retained_minus_excluded = domain_difference_ci$estimate,
    domain_difference_ci_lower = domain_difference_ci$lower,
    domain_difference_ci_upper = domain_difference_ci$upper,
    domain_difference_excludes_zero = domain_difference_ci$excludes_zero,
    n_domain_difference_blocks = as.integer(domain_difference_ci$n_blocks)
  )
})
close(pb)

lc_tab <- left_join(lc_tab, ci_list, by = "lc_id")

# write full table
out_csv_full <- file.path(
  outdir_tbl,
  sprintf("lc_fraction_class_%s_%s_%s_%s.csv", var, metric, mask, suffix)
)

write_csv(lc_tab, sub("\\.csv$", "_full_precision.csv", out_csv_full))
write_csv(round_numeric(lc_tab, 5), out_csv_full)
# paper table

paper_tab <- lc_tab |>
  transmute(
    lc_id,
    lc_name,
    area_unmasked_mkm2 = area_unm_mkm2,
    area_masked_mkm2 = area_msk_mkm2,
    area_masked_out_mkm2 = area_out_mkm2,
    retained_pct = 100 * frac_retained,
    trend_unmasked = mean_unmasked,
    trend_unmasked_ci_lower = ci_unm_lower,
    trend_unmasked_ci_upper = ci_unm_upper,
    trend_masked = mean_masked,
    trend_masked_ci_lower = ci_msk_lower,
    trend_masked_ci_upper = ci_msk_upper,
    trend_masked_out = mean_masked_out,
    trend_masked_out_ci_lower = ci_out_lower,
    trend_masked_out_ci_upper = ci_out_upper,
    retained_minus_excluded,
    difference_ci_lower,
    difference_ci_upper,
    difference_excludes_zero,
    n_difference_blocks,
    domain_retained_minus_excluded,
    domain_difference_ci_lower,
    domain_difference_ci_upper,
    domain_difference_excludes_zero,
    n_domain_difference_blocks,
    trend_delta = mean_masked - mean_unmasked,
    trend_removed = ifelse(
      (1 - frac_retained) > 0.01,
      (mean_unmasked - (frac_retained * mean_masked)) / (1 - frac_retained),
      NA_real_
    ),
    trend_delta_pct = 100 * safe_division(mean_masked - mean_unmasked, abs(mean_unmasked)),
    trend_retained_ratio = 100 * safe_division(mean_masked, mean_unmasked),
    trend_removed_pct = 100 * safe_division((
      mean_unmasked - (frac_retained * mean_masked)
    ), (1 - frac_retained)),
    removed_contrib_pct = 100 * safe_division(mean_unmasked - mean_masked, mean_unmasked)
  )

out_csv_paper <- file.path(
  outdir_tbl,
  sprintf(
    "table_lc_fraction_trend_summary_%s_%s_%s_%s_%s_main.csv",
    var,
    metric,
    mask,
    alpha,
    suffix
  )
)
write_csv(paper_tab, sub("\\.csv$", "_full_precision.csv", out_csv_paper))
write_csv(round_numeric(paper_tab, 5), out_csv_paper)
#  plot
plot_tab <- lc_tab |>
  filter(!is.na(lc_id), area_unm_mkm2 > 0, !lc_id %in% nonveg_lc_ids) |>
  mutate(
    trend_delta = mean_masked - mean_unmasked,
    frac_retained_label = sprintf("%.0f%%", frac_retained * 100)
  ) |>
  arrange(mean_masked) |>
  mutate(
    lc_name = factor(lc_name, levels = lc_name),
    sig_unmasked = ci_unm_lower * ci_unm_upper > 0,
    sig_masked = ci_msk_lower * ci_msk_upper > 0
  )
plot_long <- plot_tab |>
  select(
    lc_name,
    mean_unmasked,
    mean_masked,
    ci_unm_lower,
    ci_unm_upper,
    ci_msk_lower,
    ci_msk_upper,
    sig_unmasked,
    sig_masked
  ) |>
  pivot_longer(
    cols = c(mean_unmasked, mean_masked),
    names_to = "scenario",
    values_to = "trend"
  ) |>
  mutate(
    ci_lower = if_else(scenario == "mean_unmasked", ci_unm_lower, ci_msk_lower),
    ci_upper = if_else(scenario == "mean_unmasked", ci_unm_upper, ci_msk_upper),
    sig = if_else(scenario == "mean_unmasked", sig_unmasked, sig_masked),
    scenario = factor(
      scenario,
      levels = c("mean_unmasked", "mean_masked"),
      labels = c("Unmasked", "Masked")
    ),
    shape_type = ifelse(sig, 16, 1)
  )
if (!use_relative) {
  plot_tab <- plot_tab |>
    mutate(across(c(mean_unmasked, mean_masked, ci_unm_lower, ci_unm_upper,
                    ci_msk_lower, ci_msk_upper), ~ 1000 * .x))
  plot_long <- plot_long |>
    mutate(across(c(trend, ci_lower, ci_upper), ~ 1000 * .x))
}
p <- plot_lc_trend(plot_tab, plot_long, plot_scale_factor, unit_label)
# Save
out_png <- file.path(
  outdir_fig,
  sprintf(
    "lc_fraction_class_trend_summary_%s_%s_%s_%s_%s_main.png",
    var,
    metric,
    mask,
    alpha,
    suffix
  )
)
out_pdf <- file.path(
  outdir_fig,
  sprintf(
    "lc_fraction_class_trend_summary_%s_%s_%s_%s_%s_main.pdf",
    var,
    metric,
    mask,
    alpha,
    suffix
  )
)

ggsave(out_png,
  p,
  width = 8.4,
  height = 6.2,
  dpi = 320
)
ggsave(out_pdf,
  p,
  width = 8.4,
  height = 6.2,
  dpi = 320
)
