# ==============================================================================
# 02_stable_lc_and_trend_summary.R
# From yearly 0.25° majority LC maps -> stable LC (mode across years),
# then area-weighted trend summaries by LC class (one zonal() call).
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(here)
  library(ggplot2)
  library(tidyr)
})

terraOptions(progress = 0, memfrac = 0.6, todisk = TRUE)

# ---- config ------------------------------------------------------------------
tau <- "tau_0.1"
mask <- "CCI"
var <- "LAI"
metric <- "yearmean"
use_relative <- FALSE

lc_years <- 1992:2020
stability_min <- NULL # e.g. 0.6; NULL disables

lc_yearly_dir <- here("analysis", "tmp", "lc025_majority_yearly")
out_dir <- here("analysis", "results", "lc_trends", tau)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

merge_id <- function(x) {
  y <- merge_to_parent[as.character(x)]
  ifelse(is.na(y), x, unname(y))
}
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
stopifnot(all(file.exists(lc_files)))

lc_stack <- rast(lc_files)
lc_stack <- align_to_ref(lc_stack, ref025, method = "near") |> mask(area)
# Merge subclasses on the 0.25° yearly maps BEFORE stability
# (classification-style mapping; non-listed IDs remain unchanged)
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
unit_label <- if (use_relative) "%/yr" else "per-yr (native units)"

# ---- (3) ONE zonal() call (fast) --------------------------------------------
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

lc_tab <- zone |>
  mutate(
    mean_unmasked   = scale_factor * (num_unm / den_unm_km2),
    mean_masked     = scale_factor * (num_msk / den_msk_km2),
    mean_masked_out = scale_factor * (num_out / den_out_km2),
    area_unm_mkm2   = den_unm_km2 / 1e6,
    area_msk_mkm2   = den_msk_km2 / 1e6,
    area_out_mkm2   = den_out_km2 / 1e6,
    frac_retained   = den_msk_km2 / den_unm_km2
  ) |>
  mutate(across(starts_with("mean_"), ~ ifelse(is.finite(.x), .x, NA_real_))) |>
  arrange(desc(area_unm_mkm2))

write.csv(
  lc_tab,
  file.path(out_dir, sprintf("lc_class_%s_%s_%s_%s.csv", var, metric, mask, suffix)),
  row.names = FALSE
)

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
plot_tab <- lc_tab |>
  filter(!is.na(lc_id)) |>
  left_join(lc_legend_merged, by = "lc_id") |>
  mutate(lc_name = coalesce(lc_name, paste0("Class ", lc_id))) |>
  arrange(desc(mean_unmasked))

plot_df <- plot_tab |>
  select(lc_name, mean_unmasked, mean_masked, frac_retained) |>
  pivot_longer(c(mean_unmasked, mean_masked), names_to = "metric", values_to = "value") |>
  mutate(
    metric  = factor(metric, levels = c("mean_unmasked", "mean_masked"), labels = c("Unmasked", "Masked")),
    lc_name = factor(lc_name, levels = rev(unique(plot_tab$lc_name)))
  )

pd <- position_dodge(width = 0.7)
yspan <- diff(range(plot_df$value, na.rm = TRUE))
yoff <- 0.03 * ifelse(is.finite(yspan) && yspan > 0, yspan, max(abs(plot_df$value), na.rm = TRUE))

area_unm_tot <- sum(lc_tab$area_unm_mkm2, na.rm = TRUE)
area_msk_tot <- sum(lc_tab$area_msk_mkm2, na.rm = TRUE)

p_bar <- ggplot(plot_df, aes(lc_name, value, fill = metric)) +
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
  labs(
    y = unit_label, fill = NULL,
    title = sprintf("Stable LC classes (%d–%d) from yearly 0.25° majority", min(lc_years), max(lc_years)),
    subtitle = sprintf(
      "Total land area: unmasked %.1f Mkm²; masked %.1f Mkm²; labels = retained area (%%)",
      area_unm_tot, area_msk_tot
    )
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", plot.margin = margin(5.5, 25, 5.5, 5.5))

ggsave(
  file.path(out_dir, sprintf("plot_lc_bar_%s_%s_%s_%s.png", var, metric, mask, suffix)),
  p_bar,
  width = 8, height = 6, dpi = 300
)

cat("Wrote:\n")
cat("  ", file.path(out_dir, sprintf("lc_class_%s_%s_%s_%s.csv", var, metric, mask, suffix)), "\n", sep = "")
cat("  ", file.path(out_dir, sprintf("plot_lc_bar_%s_%s_%s_%s.png", var, metric, mask, suffix)), "\n", sep = "")
cat("  ", lc_cache, "\n", sep = "")
