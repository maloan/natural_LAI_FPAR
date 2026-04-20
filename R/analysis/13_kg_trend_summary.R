# ==============================================================================
# 13_kg_trend_summary.R
# KG summaries of trends:
#   (1) FULL KG classes (3-letter codes) written to CSV
#   (2) Plot aggregated by higher-order KG group (2-letter)
# + masked-out-only summary
# + bar plot (paper style) with retained-area (%) on masked bars
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(here)
  library(kgc)
  library(ggplot2)
  library(tidyr)
  library(readr)
})
source(here::here("R", "helpers", "cli_args.R"))
source(here::here("R", "helpers", "analysis_plot_helpers.R"))

# ---- config ------------------------------------------------------------------
default_cfg <- list(
  tau = "tau_0.1",
  mask = "CCI",
  var = "LAI",
  metric = "yearmean",
  use_relative = TRUE,
  kg_res = "course",
  chunk_size = 50000L,
  top_n = Inf,
  top_n3 = Inf,
  make_preview_map = FALSE
)

cfg <- parse_cli_args(default_cfg)

tau <- as.character(cfg$tau)
mask <- as.character(cfg$mask)
var <- as.character(cfg$var)
metric <- as.character(cfg$metric)
use_relative <- isTRUE(cfg$use_relative)
kg_res <- as.character(cfg$kg_res)
chunk_size <- as.integer(cfg$chunk_size)
top_n <- cfg$top_n
top_n3 <- cfg$top_n3
make_preview_map <- isTRUE(cfg$make_preview_map)

if (is.na(chunk_size) || chunk_size <= 0) {
  stop("chunk_size must be a positive integer")
}


# ---- inputs ------------------------------------------------------------------
ref025 <- terra::rast(here::here("src", "ref_0p25.nc"))
area <- terra::rast(here::here("src", "area_0p25_validdomain_km2.nc"))[[1]]

if (use_relative) {
  f_unm <- here::here(
    "analysis", "unmasked", "0p25",
    sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
  f_msk <- here::here(
    "output", tau, "eval", sprintf("trend_%s_%s", var, mask),
    sprintf("%s_%s_trend_relative_peryear_0p25.nc", var, metric)
  )
} else {
  f_unm <- here::here(
    "analysis", "unmasked", "0p25",
    sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, metric)
  )
  f_msk <- here::here(
    "output", tau, "eval", sprintf("trend_%s_%s", var, mask),
    sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, metric)
  )
}

out_dir <- here::here("analysis", "results", "kg_trends", tau)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
paper_fig_dir <- here::here("analysis", "results", "paper_figures")
dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)

align_to_ref <- function(r, ref, method = "near") {
  if (compareGeom(r, ref, stopOnError = FALSE)) {
    return(r)
  }
  resample(r, ref, method = method)
}

# ---- read rasters ------------------------------------------------------------
stopifnot(file.exists(f_unm), file.exists(f_msk))
r_unm <- terra::rast(f_unm)[[1]] |> align_to_ref(ref025, method = "bilinear")
r_msk <- terra::rast(f_msk)[[1]] |> align_to_ref(ref025, method = "bilinear")

# ---- units -------------------------------------------------------------------
scale_factor <- if (use_relative) 100 else 1
suffix <- if (use_relative) "rel" else "abs"
unit_label <- if (use_relative) "% yr-1" else sprintf("%s yr-1", var)

# ---- helpers -----------------------------------------------------------------
zonal_weighted_mean <- function(num, w, z, zone_name = "zone") {
  zn <- terra::zonal(num, z, fun = "sum", na.rm = TRUE)
  zw <- terra::zonal(w, z, fun = "sum", na.rm = TRUE)
  colnames(zn) <- c(zone_name, "num")
  colnames(zw) <- c(zone_name, "den_km2")
  as.data.frame(zn) |>
    dplyr::full_join(as.data.frame(zw), by = zone_name) |>
    dplyr::mutate(dplyr::across(c("num", "den_km2"), as.numeric))
}

lookup_cz_chunked <- function(pts, chunk_size = 50000L, res = "course") {
  n <- nrow(pts)
  out <- character(n)
  starts <- seq.int(1L, n, by = chunk_size)
  for (k in seq_along(starts)) {
    i1 <- starts[k]
    i2 <- min(n, i1 + chunk_size - 1L)
    cat(sprintf("LookupCZ chunk %d/%d (%d-%d)\n", k, length(starts), i1, i2))
    out[i1:i2] <- as.character(kgc::LookupCZ(pts[i1:i2, ], res = res, rc = TRUE))
  }
  out
}

kg_code2_name <- function(code2) {
  dplyr::case_when(
    code2 == "Af" ~ "Tropical rainforest",
    code2 == "Am" ~ "Tropical monsoon",
    code2 == "As" ~ "Tropical savanna (dry summer)",
    code2 == "Aw" ~ "Tropical savanna (dry winter)",
    code2 == "BW" ~ "Desert",
    code2 == "BS" ~ "Steppe",
    code2 == "Cf" ~ "Temperate, fully humid",
    code2 == "Cl" ~ "Unknown climate zone",
    code2 == "Cs" ~ "Temperate, dry summer",
    code2 == "Cw" ~ "Temperate, dry winter",
    code2 == "Df" ~ "Cold, fully humid",
    code2 == "Ds" ~ "Cold, dry summer",
    code2 == "Dw" ~ "Cold, dry winter",
    code2 == "ET" ~ "Tundra",
    code2 == "EF" ~ "Ice cap",
    TRUE ~ code2
  )
}

# ---- KG codes on grid (land only) --------------------------------------------
kg_cache <- file.path(out_dir, sprintf("kg_code_grid_%s.rds", kg_res))

if (file.exists(kg_cache)) {
  kg_code <- readRDS(kg_cache)
} else {
  valid_cells <- which(!is.na(terra::values(area)))
  xy <- terra::xyFromCell(ref025, valid_cells)
  pts <- data.frame(Site = seq_len(nrow(xy)), Longitude = xy[, 1], Latitude = xy[, 2])
  kg_code_valid <- lookup_cz_chunked(pts, chunk_size = chunk_size, res = kg_res)
  kg_code <- rep(NA_character_, terra::ncell(ref025))
  kg_code[valid_cells] <- kg_code_valid
  saveRDS(kg_code, kg_cache)
}

codes <- sort(unique(kg_code))
codes <- codes[!is.na(codes) & codes != ""]
stopifnot(length(codes) > 0)

kg_id <- ref025[[1]]
terra::values(kg_id) <- match(kg_code, codes)
names(kg_id) <- "kg_id"

# ---- weights/numerators ------------------------------------------------------
w_unm <- area
w_msk <- ifel(is.na(r_msk), NA_real_, area)
w_out <- ifel(!is.na(r_unm) & is.na(r_msk), area, NA_real_)

num_unm <- r_unm * w_unm
num_msk <- r_msk * w_msk
num_out <- r_unm * w_out

# ---- zonal sums by full code (kg_id) -----------------------------------------
z_unm <- zonal_weighted_mean(num_unm, w_unm, kg_id, "kg_id") |>
  dplyr::rename(num_unm = num, den_unm_km2 = den_km2)
z_msk <- zonal_weighted_mean(num_msk, w_msk, kg_id, "kg_id") |>
  dplyr::rename(num_msk = num, den_msk_km2 = den_km2)
z_out <- zonal_weighted_mean(num_out, w_out, kg_id, "kg_id") |>
  dplyr::rename(num_out = num, den_out_km2 = den_km2)

# ------------------------------------------------------------------------------
# (1) FULL KG classes table (3-letter codes)  -> CSV
# ------------------------------------------------------------------------------
kg_base <- z_unm |>
  dplyr::full_join(z_msk, by = "kg_id") |>
  dplyr::full_join(z_out, by = "kg_id") |>
  dplyr::mutate(
    kg_code  = codes[kg_id],
    kg_code2 = substr(kg_code, 1, 2),
    kg_code3 = substr(kg_code, 1, 3)
  )

kg_full <- kg_base |>
  dplyr::mutate(
    mean_unmasked   = (num_unm / den_unm_km2),
    mean_masked     = (num_msk / den_msk_km2),
    mean_masked_out = (num_out / den_out_km2),
    area_unm_mkm2   = den_unm_km2 / 1e6,
    area_msk_mkm2   = den_msk_km2 / 1e6,
    area_out_mkm2   = den_out_km2 / 1e6,
    frac_retained   = den_msk_km2 / den_unm_km2
  ) |>
  dplyr::mutate(
    mean_unmasked   = scale_factor * mean_unmasked,
    mean_masked     = scale_factor * mean_masked,
    mean_masked_out = scale_factor * mean_masked_out
  ) |>
  dplyr::mutate(dplyr::across(starts_with("mean_"), ~ ifelse(is.finite(.x), .x, NA_real_))) |>
  dplyr::select(
    kg_id, kg_code, kg_code2, kg_code3,
    mean_unmasked, mean_masked, mean_masked_out,
    area_unm_mkm2, area_msk_mkm2, area_out_mkm2,
    frac_retained
  ) |>
  dplyr::arrange(dplyr::desc(area_unm_mkm2))

out_csv_full <- file.path(out_dir, sprintf("kg_class_full_%s_%s_%s_%s.csv", var, metric, mask, suffix))
write_csv(kg_full, out_csv_full)

# ------------------------------------------------------------------------------
# (1b) 3-LETTER KG classes table (e.g., BWh, Cfa) -> CSV
# ------------------------------------------------------------------------------
kg_tab3 <- kg_base |>
  dplyr::group_by(kg_code3) |>
  dplyr::summarise(
    num_unm = sum(num_unm, na.rm = TRUE),
    den_unm_km2 = sum(den_unm_km2, na.rm = TRUE),
    num_msk = sum(num_msk, na.rm = TRUE),
    den_msk_km2 = sum(den_msk_km2, na.rm = TRUE),
    num_out = sum(num_out, na.rm = TRUE),
    den_out_km2 = sum(den_out_km2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    mean_unmasked   = num_unm / den_unm_km2,
    mean_masked     = num_msk / den_msk_km2,
    mean_masked_out = num_out / den_out_km2,
    area_unm_mkm2   = den_unm_km2 / 1e6,
    area_msk_mkm2   = den_msk_km2 / 1e6,
    area_out_mkm2   = den_out_km2 / 1e6,
    frac_retained   = ifelse(den_unm_km2 > 0, den_msk_km2 / den_unm_km2, NA_real_)
  ) |>
  dplyr::mutate(
    mean_unmasked   = scale_factor * mean_unmasked,
    mean_masked     = scale_factor * mean_masked,
    mean_masked_out = scale_factor * mean_masked_out
  ) |>
  dplyr::select(
    kg_code3,
    mean_unmasked, mean_masked, mean_masked_out,
    area_unm_mkm2, area_msk_mkm2, area_out_mkm2,
    frac_retained
  ) |>
  dplyr::arrange(dplyr::desc(area_unm_mkm2))

out_csv_code3 <- file.path(out_dir, sprintf("kg_code3_%s_%s_%s_%s.csv", var, metric, mask, suffix))
write_csv(kg_tab3, out_csv_code3)

# ---- bar plot (3-letter classes) ---------------------------------------------
show_n3 <- if (is.infinite(top_n3)) nrow(kg_tab3) else min(top_n3, nrow(kg_tab3))

plot_tab3 <- kg_tab3 |>
  dplyr::arrange(dplyr::desc(mean_unmasked)) |>
  dplyr::slice_head(n = show_n3) |>
  dplyr::mutate(kg_label = kg_code3)

area3_unm_tot <- sum(kg_tab3$area_unm_mkm2, na.rm = TRUE)
area3_msk_tot <- sum(kg_tab3$area_msk_mkm2, na.rm = TRUE)

plot_df3 <- plot_tab3 |>
  dplyr::select(kg_label, mean_unmasked, mean_masked, frac_retained) |>
  tidyr::pivot_longer(c(mean_unmasked, mean_masked),
    names_to = "metric", values_to = "value"
  ) |>
  dplyr::mutate(
    metric = factor(metric,
      levels = c("mean_unmasked", "mean_masked"),
      labels = c("Unmasked", "Masked")
    ),
    kg_label = factor(kg_label, levels = rev(unique(plot_tab3$kg_label)))
  )

pd3 <- position_dodge(width = 0.7)
finite_vals3 <- plot_df3$value[is.finite(plot_df3$value)]
if (!length(finite_vals3)) {
  yoff3 <- 0
} else {
  yspan3 <- diff(range(finite_vals3))
  yoff3 <- 0.03 * ifelse(is.finite(yspan3) && yspan3 > 0, yspan3, max(abs(finite_vals3)))
}

labs3 <- labs(
  y = unit_label, fill = NULL,
  title = "Köppen–Geiger classes (3-letter)",
  subtitle = sprintf(
    "Total land area: unmasked %.1f Mkm²; masked %.1f Mkm²; labels = retained area (%%)",
    area3_unm_tot, area3_msk_tot
  )
)

if (!nrow(plot_df3) || !length(finite_vals3)) {
  p_bar3 <- ggplot() +
    annotate("text", x = 1, y = 1, label = "No finite KG class trends available") +
    xlim(0.5, 1.5) +
    ylim(0.5, 1.5) +
    theme_map() +
    labs3 +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title.position = "plot"
    )
} else {
  p_bar3 <- ggplot(plot_df3, aes(x = kg_label, y = value, fill = metric)) +
    geom_col(position = pd3, width = 0.65, na.rm = TRUE) +
    geom_text(
      data = dplyr::filter(plot_df3, metric == "Masked"),
      aes(
        y = value + sign(value) * yoff3,
        label = sprintf("%.0f%%", 100 * frac_retained)
      ),
      position = pd3,
      size = 2.6,
      na.rm = TRUE
    ) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    scale_fill_manual(values = c("Unmasked" = "grey40", "Masked" = "black")) +
    labs3 +
    theme_map() +
    theme(
      axis.title.x = element_text(size = 9),
      legend.position = "top",
      plot.title.position = "plot",
      plot.margin = margin(5.5, 25, 5.5, 5.5)
    )
}

ggsave(
  file.path(out_dir, sprintf(
    "plot_kg_code3_bar_top%d_%s_%s_%s_%s.png",
    show_n3, var, metric, mask, suffix
  )),
  p_bar3,
  width = 8.0, height = 6.0, dpi = 300
)

ggsave(
  file.path(out_dir, sprintf(
    "plot_kg_code3_bar_top%d_%s_%s_%s_%s.pdf",
    show_n3, var, metric, mask, suffix
  )),
  p_bar3,
  width = 8.0, height = 6.0
)

# ------------------------------------------------------------------------------
# (2) HIGHER-ORDER plot table (2-letter groups) -> CSV + plot
# ------------------------------------------------------------------------------
kg_tab2 <- kg_base |>
  dplyr::group_by(kg_code2) |>
  dplyr::summarise(
    num_unm = sum(num_unm, na.rm = TRUE),
    den_unm_km2 = sum(den_unm_km2, na.rm = TRUE),
    num_msk = sum(num_msk, na.rm = TRUE),
    den_msk_km2 = sum(den_msk_km2, na.rm = TRUE),
    num_out = sum(num_out, na.rm = TRUE),
    den_out_km2 = sum(den_out_km2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    mean_unmasked   = num_unm / den_unm_km2,
    mean_masked     = num_msk / den_msk_km2,
    mean_masked_out = num_out / den_out_km2,
    area_unm_mkm2   = den_unm_km2 / 1e6,
    area_msk_mkm2   = den_msk_km2 / 1e6,
    area_out_mkm2   = den_out_km2 / 1e6,
    frac_retained   = ifelse(den_unm_km2 > 0, den_msk_km2 / den_unm_km2, NA_real_),
    kg_group_name   = kg_code2_name(kg_code2)
  ) |>
  dplyr::mutate(
    mean_unmasked   = scale_factor * mean_unmasked,
    mean_masked     = scale_factor * mean_masked,
    mean_masked_out = scale_factor * mean_masked_out
  ) |>
  dplyr::select(
    kg_code2, kg_group_name,
    mean_unmasked, mean_masked, mean_masked_out,
    area_unm_mkm2, area_msk_mkm2, area_out_mkm2,
    frac_retained
  ) |>
  dplyr::arrange(dplyr::desc(area_unm_mkm2))

out_csv_code2 <- file.path(out_dir, sprintf("kg_code2_%s_%s_%s_%s.csv", var, metric, mask, suffix))
write_csv(kg_tab2, out_csv_code2)

# ---- bar plot (2-letter groups) ----------------------------------------------
show_n <- if (is.infinite(top_n)) nrow(kg_tab2) else min(top_n, nrow(kg_tab2))

plot_tab <- kg_tab2 |>
  dplyr::arrange(dplyr::desc(mean_unmasked)) |>
  dplyr::slice_head(n = show_n) |>
  dplyr::mutate(kg_label = paste0(kg_code2, " — ", kg_group_name))

area_unm_tot <- sum(kg_tab2$area_unm_mkm2, na.rm = TRUE)
area_msk_tot <- sum(kg_tab2$area_msk_mkm2, na.rm = TRUE)

plot_df <- plot_tab |>
  dplyr::select(kg_label, mean_unmasked, mean_masked, frac_retained) |>
  tidyr::pivot_longer(c(mean_unmasked, mean_masked),
    names_to = "metric", values_to = "value"
  ) |>
  dplyr::mutate(
    metric = factor(metric,
      levels = c("mean_unmasked", "mean_masked"),
      labels = c("Unmasked", "Masked")
    ),
    kg_label = factor(kg_label, levels = rev(unique(plot_tab$kg_label)))
  )

pd <- position_dodge(width = 0.7)
finite_vals <- plot_df$value[is.finite(plot_df$value)]
if (!length(finite_vals)) {
  yoff <- 0
} else {
  yspan <- diff(range(finite_vals))
  yoff <- 0.03 * ifelse(is.finite(yspan) && yspan > 0, yspan, max(abs(finite_vals)))
}

labs2 <- labs(
  y = unit_label, fill = NULL,
  title = "Köppen–Geiger groups (2-letter)",
  subtitle = sprintf(
    "Total land area: unmasked %.1f Mkm²; masked %.1f Mkm²; labels = retained area (%%)",
    area_unm_tot, area_msk_tot
  )
)

if (!nrow(plot_df) || !length(finite_vals)) {
  p_bar <- ggplot() +
    annotate("text", x = 1, y = 1, label = "No finite KG group trends available") +
    xlim(0.5, 1.5) +
    ylim(0.5, 1.5) +
    theme_map() +
    labs2 +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title.position = "plot"
    )
} else {
  p_bar <- ggplot(plot_df, aes(x = kg_label, y = value, fill = metric)) +
    geom_col(position = pd, width = 0.65, na.rm = TRUE) +
    geom_text(
      data = dplyr::filter(plot_df, metric == "Masked"),
      aes(
        y = value + sign(value) * yoff,
        label = sprintf("%.0f%%", 100 * frac_retained)
      ),
      position = pd,
      size = 2.6,
      na.rm = TRUE
    ) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    scale_fill_manual(values = c("Unmasked" = "grey40", "Masked" = "black")) +
    labs2 +
    theme_map() +
    theme(
      axis.title.x = element_text(size = 9),
      legend.position = "top",
      plot.title.position = "plot",
      plot.margin = margin(5.5, 25, 5.5, 5.5)
    )
}

ggsave(
  file.path(out_dir, sprintf(
    "plot_kg_code2_bar_top%d_%s_%s_%s_%s.png",
    show_n, var, metric, mask, suffix
  )),
  p_bar,
  width = 8.0, height = 6.0, dpi = 300
)

out_pdf_code2 <- file.path(out_dir, sprintf(
  "plot_kg_code2_bar_top%d_%s_%s_%s_%s.pdf",
  show_n, var, metric, mask, suffix
))
ggsave(out_pdf_code2, p_bar, width = 8.0, height = 6.0)

paper_png <- file.path(
  paper_fig_dir,
  sprintf("kg_group_trend_summary_%s_%s_%s_%s_main.png", var, metric, mask, tau)
)
paper_pdf <- file.path(
  paper_fig_dir,
  sprintf("kg_group_trend_summary_%s_%s_%s_%s_main.pdf", var, metric, mask, tau)
)
ggsave(paper_png, p_bar, width = 8.2, height = 6.1, dpi = 320)
ggsave(paper_pdf, p_bar, width = 8.2, height = 6.1)

paper_tab <- kg_tab2 |>
  dplyr::transmute(
    kg_code2,
    kg_group_name,
    area_unmasked_mkm2 = area_unm_mkm2,
    area_masked_mkm2 = area_msk_mkm2,
    retained_pct = 100 * frac_retained,
    trend_unmasked = mean_unmasked,
    trend_masked = mean_masked,
    trend_delta = mean_masked - mean_unmasked
  ) |>
  dplyr::mutate(dplyr::across(where(is.double), ~ round(.x, 4))) |>
  dplyr::arrange(dplyr::desc(area_unmasked_mkm2))

out_csv_paper <- file.path(
  paper_fig_dir,
  sprintf("table_kg_group_trend_summary_%s_%s_%s_%s_main.csv", var, metric, mask, tau)
)
write_csv(paper_tab, out_csv_paper)

# ---- optional map of one group ----------------------------------------------
if (isTRUE(make_preview_map)) {
  target <- "Dsb"
  r_target <- ref025[[1]]
  terra::values(r_target) <- substr(kg_code, 1, 3) == target
  names(r_target) <- paste0("is_", target)
  plot(r_target, main = paste0("KG group = ", target))
}

cat("Done.\n")
cat("Wrote:\n")
cat("  ", out_csv_full, "\n", sep = "")
cat("  ", out_csv_code3, "\n", sep = "")
cat("  ", out_csv_code2, "\n", sep = "")
cat("  ", out_csv_paper, "\n", sep = "")
cat("  ", paper_png, "\n", sep = "")
cat("  ", paper_pdf, "\n", sep = "")
