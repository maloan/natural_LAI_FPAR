# ==============================================================================
# 12_a_landcover_trend_summary.R — Fractional land-cover class trend summaries
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
  progress = 1,
  memfrac = 0.6,
  todisk = TRUE
)

# config
default_cfg <- list(
  tau = "tau_0.1",
  mask = "CCI",
  var = "LAI",
  metric = "yearmean",
  use_relative = TRUE,
  lc_year_start = 1992L,
  lc_year_end = 2021L,
  n_boot = 1000L,
  conf = 0.95
)

cfg <- parse_cli_args(default_cfg)

tau <- as.character(cfg$tau)
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
lc_frac_dir <- here("analysis", "tmp", "lc025_fraction_yearly")
outdir_fig <- here("analysis", "results", "figures", "summaries")
outdir_tbl <- here("analysis", "results", "tables", "land_cover")

dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tbl, recursive = TRUE, showWarnings = FALSE)

ref025 <- terra::rast(here::here("src", "ref_0p25.nc"))
area <- rast(here::here("src", "area_0p25_validdomain_km2.nc"))[[1]]
area_vals <- terra::values(area, dataframe = FALSE)
valid_domain_cells <- which(is.finite(area_vals) & area_vals > 0)
block_id <- make_block_id(area, block_size_deg = 5)

trend_files <- function(use_relative) {
  suf <- if (use_relative) {
    "trend_relative_peryear"
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
      tau,
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

#  read and average fractional land cover
lc_files <- file.path(lc_frac_dir, sprintf("lc025_fraction_%d.tif", lc_years))
missing_lc <- lc_files[!file.exists(lc_files)]
if (length(missing_lc)) {
  stop("Missing yearly fractional LC files")
}

message("Reading fractional land-cover files...")
frac_stack <- rast(lc_files)
frac_stack <- align_to_template(frac_stack, ref025, method = "bilinear")
frac_stack <- mask(frac_stack, area)
class_names <- unique(names(frac_stack))
frac_mean_file <- file.path(outdir_fig, sprintf(
  "lc025_fraction_mean_%d-%d.tif",
  min(lc_years),
  max(lc_years)
))
frac_mean_layers <- lapply(class_names, function(nm) {
  idx <- which(names(frac_stack) == nm)

  app(frac_stack[[idx]], mean, na.rm = TRUE)
})
frac_mean <- rast(frac_mean_layers)
names(frac_mean) <- class_names
writeRaster(
  frac_mean,
  filename = frac_mean_file,
  overwrite = TRUE,
  gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"),
  datatype = "FLT4S"
)

# read trends
tf <- trend_files(use_relative)
r_unm <- rast(tf$unm)[[1]]
r_msk <- rast(tf$msk)[[1]]
r_unm <- align_to_template(r_unm, ref025, method = "bilinear")
r_msk <- align_to_template(r_msk, ref025, method = "bilinear")

if (use_relative) {
  scale_factor <- 100
  suffix <- "rel"
  unit_label <- "% yr-1"
} else {
  scale_factor <- 1
  suffix <- "abs"
  unit_label <- sprintf("%s yr-1", var)
}

#  masks and base weights
w_unm_base <- ifel(!is.na(r_unm), area, NA_real_)
w_msk_base <- ifel(!is.na(r_msk), area, NA_real_)
w_out_base <- ifel(!is.na(r_unm) & is.na(r_msk), area, NA_real_)

#  fraction-weighted class summaries
summarise_fraction_class <- function(frac, lc_id) {
  w_unm_cls <- frac * w_unm_base
  w_msk_cls <- frac * w_msk_base
  w_out_cls <- frac * w_out_base

  den_unm <- global(w_unm_cls, "sum", na.rm = TRUE)[1, 1]
  den_msk <- global(w_msk_cls, "sum", na.rm = TRUE)[1, 1]
  den_out <- global(w_out_cls, "sum", na.rm = TRUE)[1, 1]

  num_unm <- global(r_unm * w_unm_cls, "sum", na.rm = TRUE)[1, 1]
  num_msk <- global(r_msk * w_msk_cls, "sum", na.rm = TRUE)[1, 1]
  num_out <- global(r_unm * w_out_cls, "sum", na.rm = TRUE)[1, 1]

  tibble::tibble(
    lc_id = lc_id,
    mean_unmasked = scale_factor * safe_division(num_unm, den_unm),
    mean_masked = scale_factor * safe_division(num_msk, den_msk),
    mean_masked_out = scale_factor * safe_division(num_out, den_out),
    area_unm_mkm2 = den_unm / 1e6,
    area_msk_mkm2 = den_msk / 1e6,
    area_out_mkm2 = den_out / 1e6,
    frac_retained = safe_division(den_msk, den_unm)
  )
}

lc_ids <- as.integer(sub("^lc_", "", names(frac_mean)))
lc_tab <- purrr::map2_dfr(seq_len(nlyr(frac_mean)), lc_ids, function(i, id) {
  summarise_fraction_class(frac_mean[[i]], id)
}) |>
  mutate(across(starts_with("mean_"), ~ ifelse(is.finite(.x), .x, NA_real_))) |>
  left_join(lc_legend, by = "lc_id") |>
  mutate(lc_name = coalesce(lc_name, paste0("Class ", lc_id))) |>
  arrange(desc(abs(mean_masked)))

# Prepare vectors for bootstrap
ref_cells <- which(is.finite(area_vals) & area_vals > 0)
block_id_sub <- block_id[ref_cells]

# Prepare vectors for bootstrap
r_unm_full <- terra::values(r_unm, dataframe = FALSE)
r_msk_full <- terra::values(r_msk, dataframe = FALSE)
w_unm_full <- terra::values(w_unm_base, dataframe = FALSE)
w_msk_full <- terra::values(w_msk_base, dataframe = FALSE)
w_out_full <- terra::values(w_out_base, dataframe = FALSE)

# Subset all base vectors to canonical domain cells ONCE
r_unm_vals <- r_unm_full[valid_domain_cells] * scale_factor
r_msk_vals <- r_msk_full[valid_domain_cells] * scale_factor
w_unm_base_vals <- w_unm_full[valid_domain_cells]
w_msk_base_vals <- w_msk_full[valid_domain_cells]
w_out_base_vals <- w_out_full[valid_domain_cells]
block_id_sub <- block_id[valid_domain_cells]

pb <- txtProgressBar(
  min = 0,
  max = length(lc_ids),
  style = 3
)

ci_list <- purrr::map2_dfr(seq_len(nlyr(frac_mean)), lc_ids, function(i, id) {
  setTxtProgressBar(pb, i)

  frac_vals_all <- terra::values(frac_mean[[i]], dataframe = FALSE)
  frac_vals <- frac_vals_all[valid_domain_cells]

  w_unm_cls <- frac_vals * w_unm_base_vals
  w_msk_cls <- frac_vals * w_msk_base_vals
  w_out_cls <- frac_vals * w_out_base_vals

  ok_unm <- is.finite(r_unm_vals) &
    is.finite(w_unm_cls) & w_unm_cls > 0 & !is.na(block_id_sub)
  ci_unm <- bootstrap_ci_global(r_unm_vals[ok_unm], w_unm_cls[ok_unm], block_id_sub[ok_unm], n_boot, conf)
  ok_msk <- is.finite(r_msk_vals) &
    is.finite(w_msk_cls) & w_msk_cls > 0 & !is.na(block_id_sub)
  ci_msk <- bootstrap_ci_global(r_msk_vals[ok_msk], w_msk_cls[ok_msk], block_id_sub[ok_msk], n_boot, conf)

  ok_out <- is.finite(r_unm_vals) &
    is.finite(w_out_cls) & w_out_cls > 0 & !is.na(block_id_sub)
  ci_out <- bootstrap_ci_global(r_unm_vals[ok_out], w_out_cls[ok_out], block_id_sub[ok_out], n_boot, conf)

  tibble::tibble(
    lc_id = id,
    ci_unm_lower = ci_unm$lower,
    ci_unm_upper = ci_unm$upper,
    ci_msk_lower = ci_msk$lower,
    ci_msk_upper = ci_msk$upper,
    ci_out_lower = ci_out$lower,
    ci_out_upper = ci_out$upper,
    n_unm = as.integer(ci_unm$n_eff),
    n_msk = as.integer(ci_msk$n_eff)
  )
})
close(pb)

lc_tab <- left_join(lc_tab, ci_list, by = "lc_id")

# write full table
out_csv_full <- file.path(
  outdir_tbl,
  sprintf("lc_fraction_class_%s_%s_%s_%s.csv", var, metric, mask, suffix)
)

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
    "table_lc_fraction_trend_summary_%s_%s_%s_%s_main.csv",
    var,
    metric,
    mask,
    tau
  )
)
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
p <- plot_lc_trend(plot_tab, plot_long, scale_factor)
# Save
out_png <- file.path(
  outdir_fig,
  sprintf(
    "lc_fraction_class_trend_summary_%s_%s_%s_%s_%s_main.png",
    var,
    metric,
    mask,
    tau,
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
    tau,
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
