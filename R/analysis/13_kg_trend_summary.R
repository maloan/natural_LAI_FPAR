# ==============================================================================
# 13_kg_trend_summary.R
# KG summaries of trends:
#   (1) FULL KG classes (3-letter codes) written to CSV
#   (2) Plot aggregated by higher-order KG group (2-letter)
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
  kg_res = "coarse",
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

paper_fig_dir <- here::here("analysis", "results", "figures", "summaries")
dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers -----------------------------------------------------------------
align_to_ref <- function(r, ref, method = "bilinear") {
  if (compareGeom(r, ref, stopOnError = FALSE)) {
    return(r)
  }

  resample(r, ref, method = method)
}

safe_divide <- function(num, den) {
  ifelse(is.finite(den) & den > 0, num / den, NA_real_)
}

zonal_weighted_mean <- function(num, w, z, zone_name = "zone") {
  zn <- terra::zonal(num, z, fun = "sum", na.rm = TRUE)
  zw <- terra::zonal(w, z, fun = "sum", na.rm = TRUE)

  colnames(zn) <- c(zone_name, "num")
  colnames(zw) <- c(zone_name, "den_km2")

  as.data.frame(zn) |>
    dplyr::full_join(as.data.frame(zw), by = zone_name) |>
    dplyr::mutate(
      dplyr::across(c("num", "den_km2"), as.numeric)
    )
}

lookup_cz_chunked <- function(pts, chunk_size = 50000L, res = "coarse") {
  n <- nrow(pts)
  out <- character(n)

  starts <- seq.int(1L, n, by = chunk_size)

  for (k in seq_along(starts)) {
    i1 <- starts[k]
    i2 <- min(n, i1 + chunk_size - 1L)

    cat(sprintf(
      "LookupCZ chunk %d/%d (%d-%d)\n",
      k, length(starts), i1, i2
    ))

    out[i1:i2] <- as.character(
      kgc::LookupCZ(
        pts[i1:i2, ],
        res = res,
        rc = TRUE
      )
    )
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

# ---- read rasters ------------------------------------------------------------
stopifnot(file.exists(f_unm), file.exists(f_msk))

r_unm <- terra::rast(f_unm)[[1]] |>
  align_to_ref(ref025, method = "bilinear")

r_msk <- terra::rast(f_msk)[[1]] |>
  align_to_ref(ref025, method = "bilinear")

# ---- units -------------------------------------------------------------------
scale_factor <- 1

suffix <- if (use_relative) {
  "rel"
} else {
  "abs"
}

unit_label <- if (use_relative) {
  "% yr-1"
} else {
  sprintf("%s yr-1", var)
}

# ---- KG codes on grid (land only) --------------------------------------------
kg_cache <- file.path(
  out_dir,
  sprintf("kg_code_grid_%s.rds", kg_res)
)

if (file.exists(kg_cache)) {
  kg_code <- readRDS(kg_cache)
} else {
  valid_cells <- which(!is.na(terra::values(area)))

  xy <- terra::xyFromCell(ref025, valid_cells)

  pts <- data.frame(
    Site = seq_len(nrow(xy)),
    Longitude = xy[, 1],
    Latitude = xy[, 2]
  )

  kg_code_valid <- lookup_cz_chunked(
    pts,
    chunk_size = chunk_size,
    res = kg_res
  )

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

# ---- zonal sums --------------------------------------------------------------
z_unm <- zonal_weighted_mean(num_unm, w_unm, kg_id, "kg_id") |>
  dplyr::rename(
    num_unm = num,
    den_unm_km2 = den_km2
  )

z_msk <- zonal_weighted_mean(num_msk, w_msk, kg_id, "kg_id") |>
  dplyr::rename(
    num_msk = num,
    den_msk_km2 = den_km2
  )

z_out <- zonal_weighted_mean(num_out, w_out, kg_id, "kg_id") |>
  dplyr::rename(
    num_out = num,
    den_out_km2 = den_km2
  )

kg_base <- z_unm |>
  dplyr::full_join(z_msk, by = "kg_id") |>
  dplyr::full_join(z_out, by = "kg_id") |>
  dplyr::mutate(
    kg_code = codes[kg_id],
    kg_code2 = substr(kg_code, 1, 2),
    kg_code3 = substr(kg_code, 1, 3)
  )

# ------------------------------------------------------------------------------
# Full KG classes
# ------------------------------------------------------------------------------
kg_full <- kg_base |>
  dplyr::mutate(
    mean_unmasked = safe_divide(num_unm, den_unm_km2),
    mean_masked = safe_divide(num_msk, den_msk_km2),
    mean_masked_out = safe_divide(num_out, den_out_km2),
    area_unm_mkm2 = den_unm_km2 / 1e6,
    area_msk_mkm2 = den_msk_km2 / 1e6,
    area_out_mkm2 = den_out_km2 / 1e6,
    frac_retained = safe_divide(den_msk_km2, den_unm_km2)
  ) |>
  dplyr::mutate(
    mean_unmasked   = scale_factor * mean_unmasked,
    mean_masked     = scale_factor * mean_masked,
    mean_masked_out = scale_factor * mean_masked_out
  ) |>
  dplyr::select(
    kg_id,
    kg_code,
    kg_code2,
    kg_code3,
    mean_unmasked,
    mean_masked,
    mean_masked_out,
    area_unm_mkm2,
    area_msk_mkm2,
    area_out_mkm2,
    frac_retained
  ) |>
  dplyr::arrange(dplyr::desc(area_unm_mkm2))
