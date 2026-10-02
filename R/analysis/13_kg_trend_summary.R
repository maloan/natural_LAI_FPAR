# ==============================================================================
# 13_kg_trend_summary.R — Köppen-Geiger climate-zone trend summaries
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(here)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
  library(readr)
})

source(here::here("R", "helpers", "cli_args.R"))
source(here::here("R", "helpers", "bootstrap_ci.R"))
source(here::here("R", "helpers", "plotting.R"))
source(here::here("R", "helpers", "io.R"))
source(here::here("R", "helpers", "kg_classification.R"))

default_cfg <- list(
  alpha = "alpha_0.1",
  mask = "CCI",
  var = "LAI",
  metric = "yearmean",
  use_relative = FALSE,
  n_boot = 1000L,
  conf = 0.95
)

run_kg_summary <- function(cfg) {
alpha <- as.character(cfg$alpha)
mask <- as.character(cfg$mask)
var <- as.character(cfg$var)
metric <- as.character(cfg$metric)
use_relative <- isTRUE(cfg$use_relative)
n_boot <- as.integer(cfg$n_boot)
conf <- as.numeric(cfg$conf)
stopifnot(n_boot >= 2L, conf > 0, conf < 1)
message("Climate summary: nominal 100-arc-second classification, relative=", use_relative)
terraOptions(progress = 0)

ref025 <- terra::rast(here::here("src", "ref_0p25.nc"))
area_land <- rast(here::here("src", "area_0p25_validdomain_km2.nc"))[[1]]
area <- load_summary_area(template = area_land)

land_vals <- terra::values(area_land, dataframe = FALSE)
area_vals <- terra::values(area, dataframe = FALSE)
valid_domain_cells <- which(is.finite(land_vals) & land_vals > 0)
block_id <- make_block_id(area, block_size_deg = 5)

if (use_relative) {
  trend_suffix <- "trend_relative_percent_peryear"
  suffix <- "rel"
  scale_factor <- 1
  unit_label <- expression("LAI trend (% yr"^-1 * ")")
} else {
  trend_suffix <- "trend_slope_peryear"
  suffix <- "abs"
  scale_factor <- 1
  unit_label <- expression("LAI trend (" ~ "×" ~ 10^-3 ~ m^2 ~ m^-2 ~ yr^-1 * ")")
}

f_unm <- here::here(
  "analysis",
  "unmasked",
  "0p25",
  sprintf("%s_georef_%s_%s_0p25.nc", var, metric, trend_suffix)
)

f_msk <- here::here(
  "output",
  alpha,
  "eval",
  sprintf("trend_%s_%s", var, mask),
  sprintf("%s_%s_%s_0p25.nc", var, metric, trend_suffix)
)

outdir_fig <- kg_output_dir("figures")
outdir_tbl <- kg_output_dir("tables")
dir_kg_grid <- here::here("analysis", "results", "tmp", "kg_grid")

dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_kg_grid, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(f_unm)) {
  stop("Missing unmasked trend file: ", f_unm, call. = FALSE)
}

if (!file.exists(f_msk)) {
  stop("Missing masked trend file: ", f_msk, call. = FALSE)
}

kg_code2_name <- local({
  kg2_names <- c(
    Af = "Tropical rainforest",
    Am = "Tropical monsoon",
    As = "Tropical savanna, dry summer",
    Aw = "Tropical savanna, dry winter",
    BW = "Desert",
    BS = "Steppe",
    Cf = "Temperate, fully humid",
    Cs = "Temperate, dry summer",
    Cw = "Temperate, dry winter",
    Df = "Cold, fully humid",
    Ds = "Cold, dry summer",
    Dw = "Cold, dry winter",
    ET = "Tundra",
    EF = "Ice cap"
  )
  function(code2) {
    out <- unname(kg2_names[code2])
    ifelse(is.na(out) | out == "NA", code2, out)
  }
})

make_ci <- function(r_values,
                    z_values,
                    w_values,
                    block_id,
                    n_boot,
                    conf,
                    role = c("unm", "msk")) {
  role <- match.arg(role)

  ci <- bootstrap_ci_by_class(
    r_values = r_values,
    z_values = z_values,
    w_values = w_values,
    block_id = block_id,
    n_boot = n_boot,
    conf = conf
  )

  if (role == "unm") {
    dplyr::rename(
      ci,
      n_unm = n_pixels,
      ci_unm_lower = ci_lower,
      ci_unm_upper = ci_upper,
      sig_unmasked = sig_flag
    )
  } else {
    dplyr::rename(
      ci,
      n_msk = n_pixels,
      ci_msk_lower = ci_lower,
      ci_msk_upper = ci_upper,
      sig_masked = sig_flag
    )
  }
}

make_kg_raster <- function(ref, kg_code, codes, name) {
  out <- ref[[1]]
  terra::values(out) <- match(kg_code, codes)
  names(out) <- name
  out
}

r_unm <- terra::rast(f_unm)[[1]]
r_msk <- terra::rast(f_msk)[[1]]

r_unm <- align_to_template(r_unm, ref025, method = "bilinear")
r_msk <- align_to_template(r_msk, ref025, method = "bilinear")
area_msk <- load_summary_area(template = area)

stopifnot(
  terra::compareGeom(ref025, r_unm, stopOnError = TRUE),
  terra::compareGeom(ref025, r_msk, stopOnError = TRUE),
  terra::compareGeom(ref025, area, stopOnError = TRUE),
  terra::compareGeom(ref025, area_msk, stopOnError = TRUE)
)

block_id_sub <- block_id[valid_domain_cells]

classification <- kg_grid(
  ref025,
  valid_domain_cells,
  resolution = "fine",
  cache_dir = dir_kg_grid
)
kg_code <- classification$code

coverage <- dplyr::bind_rows(lapply(c("unmasked", "masked"), function(role) {
  x <- if (role == "unmasked") values(r_unm, mat = FALSE) else values(r_msk, mat = FALSE)
  w <- if (role == "unmasked") area_vals else values(area_msk, mat = FALSE)
  valid <- is.finite(x) & is.finite(w) & w > 0
  classified <- valid & !is.na(kg_code)
  tibble(role = role, valid_cells = sum(valid), classified_cells = sum(classified),
         valid_area_km2 = sum(w[valid]), classified_area_km2 = sum(w[classified]),
         unclassified_area_pct = 100 * sum(w[valid & !classified]) / sum(w[valid]),
         classified_mean = weighted.mean(x[classified], w[classified]))
}))
write_csv(coverage, file.path(outdir_tbl, paste0("kg_coverage_", var, "_", metric, "_", suffix, ".csv")))

codes3 <- sort(unique(kg_code))
codes3 <- codes3[!is.na(codes3) & codes3 != ""]

if (!length(codes3)) {
  stop("No valid Köppen-Geiger classes found.", call. = FALSE)
}

kg_code2 <- substr(kg_code, 1, 2)
codes2 <- sort(unique(kg_code2))
codes2 <- codes2[!is.na(codes2) & codes2 != ""]

kg3_id <- make_kg_raster(ref025, kg_code, codes3, "kg3_id")
kg2_id <- make_kg_raster(ref025, kg_code2, codes2, "kg2_id")
# Valid coarse-cell estimates represent the fixed post-nonvegetated support.
w_unm <- terra::ifel(!is.na(r_unm), area, NA_real_)
w_msk <- terra::ifel(!is.na(r_msk), area_msk, NA_real_)
# The excluded-domain trend is defined only for fully excluded coarse cells.
# Excluded fine-cell trends within partially retained cells are unavailable.
w_out <- terra::ifel(!is.na(r_unm) & is.na(r_msk), area, NA_real_)

num_unm <- r_unm * w_unm
num_msk <- r_msk * w_msk
num_out <- r_unm * w_out

# Extract values once for both zone resolutions.
vals_r_unm <- terra::values(r_unm, dataframe = FALSE)
vals_r_msk <- terra::values(r_msk, dataframe = FALSE)
vals_w_unm <- terra::values(w_unm, dataframe = FALSE)
vals_w_msk <- terra::values(w_msk, dataframe = FALSE)
vals_kg3 <- terra::values(kg3_id, dataframe = FALSE)
vals_kg2 <- terra::values(kg2_id, dataframe = FALSE)

r_unm_vals <- vals_r_unm[valid_domain_cells] * scale_factor
r_msk_vals <- vals_r_msk[valid_domain_cells] * scale_factor
w_unm_vals <- vals_w_unm[valid_domain_cells]
w_msk_vals <- vals_w_msk[valid_domain_cells]
kg3_vals <- vals_kg3[valid_domain_cells]
kg2_vals <- vals_kg2[valid_domain_cells]

# Run bootstraps
ci3_unm <- make_ci(r_unm_vals,
  kg3_vals,
  w_unm_vals,
  block_id_sub,
  n_boot,
  conf,
  role = "unm"
)
ci3_msk <- make_ci(r_msk_vals,
  kg3_vals,
  w_msk_vals,
  block_id_sub,
  n_boot,
  conf,
  role = "msk"
)
ci2_unm <- make_ci(r_unm_vals,
  kg2_vals,
  w_unm_vals,
  block_id_sub,
  n_boot,
  conf,
  role = "unm"
)
ci2_msk <- make_ci(r_msk_vals,
  kg2_vals,
  w_msk_vals,
  block_id_sub,
  n_boot,
  conf,
  role = "msk"
)
kg3_codes <- tibble(
  kg3_id = seq_along(codes3),
  kg_code = codes3,
  kg_code2 = substr(codes3, 1, 2)
)

kg2_codes <- tibble(
  kg2_id = seq_along(codes2),
  kg_code2 = codes2,
  kg_name = kg_code2_name(codes2)
)
summarise_zone <- function(zone_raster,
                           code_table,
                           zone_col,
                           ci_unm,
                           ci_msk) {
  zone_name <- rlang::as_string(rlang::ensym(zone_col))
  # zonal aggregation
  out <- terra::zonal(
    c(num_unm, num_msk, num_out, w_unm, w_msk, w_out),
    zone_raster,
    fun = "sum",
    na.rm = TRUE
  )
  out <- as.data.frame(out)
  colnames(out) <- c(
    zone_name,
    "num_unm",
    "num_msk",
    "num_out",
    "area_unm_km2",
    "area_msk_km2",
    "area_fully_excluded_km2"
  )
  stopifnot(zone_name %in% names(out))
  stopifnot(all(c("num_unm", "area_unm_km2") %in% names(out)))
  # join zone metadata
  out <- left_join(out, code_table, by = zone_name)
  ci_unm_join <- ci_unm |>
    dplyr::select(-mean_est) |>
    dplyr::rename(!!zone_name := class) |>
    dplyr::mutate(!!zone_name := as.integer(.data[[zone_name]]))
  ci_msk_join <- ci_msk |>
    dplyr::select(-mean_est) |>
    dplyr::rename(!!zone_name := class) |>
    dplyr::mutate(!!zone_name := as.integer(.data[[zone_name]]))
  out <- dplyr::left_join(out, ci_unm_join, by = zone_name)
  out <- dplyr::left_join(out, ci_msk_join, by = zone_name)
  dplyr::mutate(
    out,
    mean_unmasked = scale_factor * safe_division(num_unm, area_unm_km2, positive_denominator = TRUE),
    mean_masked = scale_factor * safe_division(num_msk, area_msk_km2, positive_denominator = TRUE),
    mean_masked_out = scale_factor * safe_division(
      num_out,
      area_fully_excluded_km2,
      positive_denominator = TRUE
    ),
    diff_masked_minus_unmasked = mean_masked - mean_unmasked,
    area_unm_mkm2 = area_unm_km2 / 1e6,
    area_msk_mkm2 = area_msk_km2 / 1e6,
    area_out_mkm2 = pmax(area_unm_km2 - area_msk_km2, 0) / 1e6,
    frac_retained = safe_division(area_msk_km2, area_unm_km2, positive_denominator = TRUE),
    area_removed_pct = 100 * (1 - frac_retained),
    trend_delta_pct = 100 * safe_division(
      mean_masked - mean_unmasked,
      abs(mean_unmasked),
      positive_denominator = TRUE
    ),
    trend_retained_ratio = 100 * safe_division(mean_masked, mean_unmasked, positive_denominator = TRUE),
    trend_removed_pct = 100 * safe_division(
      mean_unmasked - (frac_retained * mean_masked),
      1 - frac_retained,
      positive_denominator = TRUE
    ),
    removed_contrib_pct = 100 * safe_division(
      mean_masked - mean_unmasked,
      mean_unmasked,
      positive_denominator = TRUE
    )
  )
}

kg_full <- summarise_zone(kg3_id, kg3_codes, kg3_id, ci3_unm, ci3_msk) |>
  select(
    kg3_id,
    kg_code,
    kg_code2,
    mean_unmasked,
    ci_unm_lower,
    ci_unm_upper,
    sig_unmasked,
    mean_masked,
    ci_msk_lower,
    ci_msk_upper,
    sig_masked,
    diff_masked_minus_unmasked,
    mean_masked_out,
    area_unm_mkm2,
    area_msk_mkm2,
    area_out_mkm2,
    area_removed_pct,
    frac_retained,
    trend_delta_pct,
    trend_retained_ratio,
    trend_removed_pct,
    removed_contrib_pct,
    n_unm,
    n_msk
  ) |>
  arrange(desc(abs(mean_masked)))

kg2_summary <- summarise_zone(kg2_id, kg2_codes, kg2_id, ci2_unm, ci2_msk) |>
  select(
    kg2_id,
    kg_code2,
    kg_name,
    mean_unmasked,
    ci_unm_lower,
    ci_unm_upper,
    sig_unmasked,
    mean_masked,
    ci_msk_lower,
    ci_msk_upper,
    sig_masked,
    diff_masked_minus_unmasked,
    mean_masked_out,
    area_unm_mkm2,
    area_msk_mkm2,
    area_out_mkm2,
    area_removed_pct,
    frac_retained,
    trend_delta_pct,
    trend_retained_ratio,
    trend_removed_pct,
    removed_contrib_pct,
    n_unm,
    n_msk
  ) |>
  arrange(desc(abs(mean_masked)))

out_full <- file.path(
  outdir_tbl,
  sprintf("kg_full_%s_%s_%s_%s.csv", var, metric, mask, suffix)
)

out_kg2 <- file.path(
  outdir_tbl,
  sprintf("kg2_summary_%s_%s_%s_%s.csv", var, metric, mask, suffix)
)

add_provenance <- function(x) mutate(x,
  kg_nominal_degrees = classification$nominal_degrees, kg_grid_signature = classification$signature,
  kgc_version = as.character(packageVersion("kgc")), bootstrap_replicates = n_boot,
  bootstrap_confidence = conf, bootstrap_block_degrees = 5)
kg_full <- add_provenance(kg_full)
kg2_summary <- add_provenance(kg2_summary)
write_csv(kg_full, out_full)
write_csv(kg2_summary, out_kg2)

valid_kg2 <- c(
  "Af",
  "Am",
  "As",
  "Aw",
  "BS",
  "BW",
  "Cf",
  "Cs",
  "Cw",
  "Df",
  "Ds",
  "Dw",
  "ET"
)

plot_tab <- kg2_summary |>
  filter(
    kg_code2 %in% valid_kg2,
    is.finite(mean_unmasked) | is.finite(mean_masked),
    is.finite(area_unm_mkm2),
    area_unm_mkm2 > 0
  ) |>
  mutate(
    kg_label = paste0(kg_code2, " - ", kg_name),
    mean_removed = abs(mean_masked - mean_unmasked),
    area_removed_label = sprintf("%.0f%%", area_removed_pct),
    frac_retained_label = sprintf("%.0f%%", frac_retained * 100),
    # Determine significance
    sig_unmasked = ci_unm_lower * ci_unm_upper > 0,
    sig_masked = ci_msk_lower * ci_msk_upper > 0
  ) |>
  arrange(desc(abs(mean_removed)))

if (!nrow(plot_tab)) {
  stop("No finite Köppen-Geiger groups available for plotting.",
    call. = FALSE
  )
}

plot_long <- plot_tab |>
  select(
    kg_label,
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

kg_levels <- rev(unique(plot_tab$kg_label))
plot_tab <- plot_tab |>
  mutate(kg_label = factor(kg_label, levels = kg_levels))
plot_long <- plot_long |>
  mutate(kg_label = factor(kg_label, levels = kg_levels))

if (!use_relative) {
  plot_tab <- plot_tab |>
    mutate(across(c(mean_unmasked, mean_masked, ci_unm_lower, ci_unm_upper,
                    ci_msk_lower, ci_msk_upper), ~ 1000 * .x))
  plot_long <- plot_long |>
    mutate(across(c(trend, ci_lower, ci_upper), ~ 1000 * .x))
}
p <- plot_kg(plot_tab, plot_long, unit_label)

out_png <- file.path(
  outdir_fig,
  sprintf(
    "kg_group_trend_summary_%s_%s_%s_%s_%s_main.png",
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
    "kg_group_trend_summary_%s_%s_%s_%s_%s_main.pdf",
    var,
    metric,
    mask,
    alpha,
    suffix
  )
)

ggsave(out_png,
  p,
  width = 8.2,
  height = 5.4,
  dpi = 320
)
ggsave(out_pdf, p, width = 8.2, height = 5.4)
invisible(list(subzone = kg_full, group = kg2_summary, coverage = coverage))
}

cfg <- parse_cli_args(default_cfg)
run_kg_summary(cfg)
