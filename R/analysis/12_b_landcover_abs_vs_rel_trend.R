# ==============================================================================
# 12_b_landcover_abs_vs_rel_trend.R — Fraction-weighted land-cover class
# absolute vs relative LAI trends
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(tibble)
})

source(here::here("R", "helpers", "cli_args.R"))
source(here::here("R", "helpers", "plotting.R"))

default_cfg <- list(
  var = "LAI",
  metric = "yearmean",
  mask = "CCI",
  tau = "tau_0.1",
  include_nonveg = FALSE,
  n_label_area = 20L,
  n_label_outlier = 5L
)

cfg <- parse_cli_args(default_cfg)

var <- as.character(cfg$var)
metric <- as.character(cfg$metric)
mask <- as.character(cfg$mask)
tau <- as.character(cfg$tau)

include_nonveg <- as.logical(cfg$include_nonveg)
n_label_area <- as.integer(cfg$n_label_area)
n_label_outlier <- as.integer(cfg$n_label_outlier)

outdir_fig <- here::here("analysis", "results", "figures", "summaries")
outdir_tbl <- here::here("analysis", "results", "tables", "land_cover")
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tbl, recursive = TRUE, showWarnings = FALSE)

summary_csv <- function(kind) {
  here::here(
    outdir_tbl,
    sprintf("lc_fraction_class_%s_%s_%s_%s.csv", var, metric, mask, kind)
  )
}

abs_csv <- summary_csv("abs")
rel_csv <- summary_csv("rel")

if (!file.exists(abs_csv)) {
  stop("Missing absolute-trend fraction summary CSV: ",
    abs_csv,
    call. = FALSE
  )
}

if (!file.exists(rel_csv)) {
  stop("Missing relative-trend fraction summary CSV: ",
    rel_csv,
    call. = FALSE
  )
}

abs_tab <- read_csv(abs_csv, show_col_types = FALSE)
rel_tab <- read_csv(rel_csv, show_col_types = FALSE)

required_abs <- c("lc_id", "mean_masked", "area_msk_mkm2")
required_rel <- c("lc_id", "mean_masked")

missing_abs <- setdiff(required_abs, names(abs_tab))
missing_rel <- setdiff(required_rel, names(rel_tab))

if (length(missing_abs) > 0) {
  stop(
    "Absolute-trend CSV is missing required columns: ",
    paste(missing_abs, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_rel) > 0) {
  stop(
    "Relative-trend CSV is missing required columns: ",
    paste(missing_rel, collapse = ", "),
    call. = FALSE
  )
}

lc_labels <- tibble::tribble(
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
# Exclude classes where LAI trends are not ecologically interpretable.
nonveg_ids <- c(0L, 190L, 200L, 210L, 220L)

plot_df <- abs_tab |>
  filter(!is.na(lc_id)) |>
  transmute(
    lc_id = as.integer(lc_id),
    trend_abs = mean_masked,
    area_mkm2 = area_msk_mkm2
  ) |>
  inner_join(
    rel_tab |>
      filter(!is.na(lc_id)) |>
      transmute(lc_id = as.integer(lc_id), trend_rel = mean_masked),
    by = "lc_id"
  ) |>
  left_join(lc_labels, by = "lc_id") |>
  mutate(lc_name = coalesce(lc_name, paste0("Class ", lc_id))) |>
  filter(
    is.finite(trend_abs),
    is.finite(trend_rel),
    is.finite(area_mkm2),
    area_mkm2 > 0
  )

if (!include_nonveg) {
  plot_df <- plot_df |>
    filter(!lc_id %in% nonveg_ids)
}

if (nrow(plot_df) == 0) {
  stop("No finite land-cover class values available for plotting.",
    call. = FALSE
  )
}

label_df <- plot_df |>
  mutate(
    area_rank = min_rank(desc(area_mkm2)),
    trend_abs_z = as.numeric(scale(trend_abs)),
    trend_rel_z = as.numeric(scale(trend_rel)),
    outlier_score = abs(trend_abs_z) + abs(trend_rel_z),
    outlier_rank = min_rank(desc(outlier_score))
  ) |>
  filter(area_rank <= n_label_area |
    outlier_rank <= n_label_outlier)

plot_df <- plot_df |>
  mutate(
    veg_type = case_when(
      lc_id %in% c(50, 60, 70, 80, 90, 160) ~ "Tree Cover",
      lc_id %in% c(100, 120) ~ "Shrubland",
      lc_id %in% c(130, 140, 150, 180) ~ "Herbaceous",
      lc_id == 10 ~ "Cropland",
      TRUE ~ "Other"
    ),
    veg_type = factor(
      veg_type,
      levels = c("Tree Cover", "Shrubland", "Herbaceous", "Cropland", "Other")
    )
  )

p <- plot_lc_abs_vs_rel(plot_df)

out_png <- file.path(
  outdir_fig,
  sprintf(
    "plot_lc_fraction_trend_abs_vs_rel_%s_%s_%s_%s.png",
    tau,
    var,
    metric,
    mask
  )
)

out_pdf <- file.path(
  outdir_fig,
  sprintf(
    "plot_lc_fraction_trend_abs_vs_rel_%s_%s_%s_%s.pdf",
    tau,
    var,
    metric,
    mask
  )
)

ggsave(out_png,
  p,
  width = 8.5,
  height = 6.1,
  dpi = 320
)
ggsave(out_pdf, p, width = 8.5, height = 6.1)
