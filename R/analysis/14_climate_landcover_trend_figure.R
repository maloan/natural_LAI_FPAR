# ==============================================================================
# 14_climate_landcover_trend_figure.R — Climate-zone and land-cover trends
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(patchwork)
  library(readr)
  library(tidyr)
})

source(here::here("R", "helpers", "cli_args.R"))
source(here::here("R", "helpers", "kg_classification.R"))
source(here::here("R", "helpers", "plotting.R"))
cfg <- parse_cli_args(list(mask = "CCI", alpha = "alpha_0.1"))
kg_dir <- kg_output_dir("tables")
lc_dir <- here::here("analysis", "results", "tables", "land_cover")
figure_dir <- kg_output_dir("figures")

kg_abs_file <- file.path(kg_dir, paste0("kg2_summary_LAI_yearmean_", cfg$mask, "_abs.csv"))
kg_rel_file <- file.path(kg_dir, paste0("kg2_summary_LAI_yearmean_", cfg$mask, "_rel.csv"))
lc_abs_file <- file.path(lc_dir, paste0("lc_fraction_class_LAI_yearmean_", cfg$mask, "_abs.csv"))
lc_rel_file <- file.path(lc_dir, paste0("lc_fraction_class_LAI_yearmean_", cfg$mask, "_rel.csv"))

input_files <- c(kg_abs_file, kg_rel_file, lc_abs_file, lc_rel_file)
if (any(!file.exists(input_files))) {
  stop("Missing input files:\n", paste(input_files[!file.exists(input_files)], collapse = "\n"))
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

kg_keep <- c("Af", "Am", "As", "Aw", "BS", "BW", "Cf", "Cs", "Cw", "Df", "Ds", "Dw", "ET")
lc_keep <- c(10L, 50L, 60L, 70L, 80L, 90L, 100L, 120L, 130L, 140L, 150L, 160L, 180L, 190L, 200L)

prepare_kg <- function(path, metric) {
  readr::read_csv(path, show_col_types = FALSE) |>
    filter(kg_code2 %in% kg_keep) |>
    mutate(kg_name = sub("^Cold", "Continental", kg_name)) |>
    transmute(
      domain = "Climate zone",
      metric,
      class = kg_code2,
      class_name = paste0(kg_code2, " - ", kg_name),
      unmasked = mean_unmasked,
      masked = mean_masked,
      unmasked_lower = ci_unm_lower,
      unmasked_upper = ci_unm_upper,
      masked_lower = ci_msk_lower,
      masked_upper = ci_msk_upper,
      retention_pct = 100 * frac_retained
    )
}

prepare_lc <- function(path, metric) {
  readr::read_csv(path, show_col_types = FALSE) |>
    filter(lc_id %in% lc_keep) |>
    transmute(
      domain = "Land-cover class",
      metric,
      class = as.character(lc_id),
      class_name = unname(c(`10` = "Cropland", `50` = "Evergreen broadleaf forest", `60` = "Deciduous broadleaf forest", `70` = "Needleleaf evergreen forest", `80` = "Needleleaf deciduous forest", `90` = "Mixed forest", `100` = "Mosaic tree/shrub", `120` = "Shrubland", `130` = "Grassland", `140` = "Lichen and mosses", `150` = "Sparse vegetation", `160` = "Flooded tree cover", `180` = "Flooded shrub/herbaceous", `190` = "Urban areas", `200` = "Bare areas")[as.character(lc_id)]),
      unmasked = mean_unmasked,
      masked = mean_masked,
      unmasked_lower = ci_unm_lower,
      unmasked_upper = ci_unm_upper,
      masked_lower = ci_msk_lower,
      masked_upper = ci_msk_upper,
      retention_pct = 100 * frac_retained
    )
}

figure_data <- bind_rows(
  prepare_kg(kg_abs_file, "Absolute"),
  prepare_kg(kg_rel_file, "Relative"),
  prepare_lc(lc_abs_file, "Absolute"),
  prepare_lc(lc_rel_file, "Relative")
) |>
  mutate(
    change_pct = 100 * (masked - unmasked) / abs(unmasked),
    retained_label = sprintf("%s  [%.0f%%]", class_name, retention_pct)
  )

figure_values_file <- here::here(
  "analysis",
  "results",
  "tables",
  "climate_landcover_trend_figure_values.csv"
)
write_csv(figure_data, figure_values_file)
write_csv(
  figure_data,
  file.path(figure_dir, "climate_landcover_trend_figure_values.csv")
)


p_kg_abs <- plot_class_trend_panel(
  filter(figure_data, domain == "Climate zone", metric == "Absolute"),
  "Absolute trend",
  expression("LAI trend (" ~ "×" ~ 10^-3 ~ m^2 ~ m^-2 ~ yr^-1 * ")"),
  "(a)",
  scale_factor = 1000
)

p_kg_rel <- plot_class_trend_panel(
  filter(figure_data, domain == "Climate zone", metric == "Relative"),
  "Relative trend",
  expression("LAI trend (% yr"^-1 * ")"),
  "(b)"
)

p_lc_abs <- plot_class_trend_panel(
  filter(figure_data, domain == "Land-cover class", metric == "Absolute"),
  "Absolute trend",
  expression("LAI trend (" ~ "×" ~ 10^-3 ~ m^2 ~ m^-2 ~ yr^-1 * ")"),
  "(a)",
  scale_factor = 1000
)

p_lc_rel <- plot_class_trend_panel(
  filter(figure_data, domain == "Land-cover class", metric == "Relative"),
  "Relative trend",
  expression("LAI trend (% yr"^-1 * ")"),
  "(b)"
)

climate_figure <- (p_kg_abs | p_kg_rel) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

landcover_figure <- (p_lc_abs | p_lc_rel) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Keep the four-panel alternative for diagnostic use.
combined_figure <- ((p_kg_abs | p_kg_rel) / ((p_lc_abs + labs(tag = "(c)")) | (p_lc_rel + labs(tag = "(d)")))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

climate_stem <- file.path(
  figure_dir,
  paste0("climate_zone_trends_absolute_relative_", cfg$mask, "_", cfg$alpha)
)
landcover_stem <- file.path(
  figure_dir,
  paste0("landcover_trends_absolute_relative_", cfg$mask, "_", cfg$alpha)
)

save_figure <- function(stem, plot, width, height) {
  ggsave(
    paste0(stem, ".png"),
    plot,
    width = width,
    height = height,
    dpi = 320
  )
  ggsave(
    paste0(stem, ".pdf"),
    plot,
    width = width,
    height = height
  )
}

save_figure(climate_stem, climate_figure, width = 13.0, height = 6.5)
save_figure(landcover_stem, landcover_figure, width = 13.0, height = 7.2)

combined_stem <- file.path(
  figure_dir,
  paste0("climate_landcover_trends_absolute_relative_", cfg$mask, "_", cfg$alpha)
)
save_figure(combined_stem, combined_figure, width = 13.0, height = 13.4)

message("Wrote separate and combined climate-zone/land-cover trend figures.")
