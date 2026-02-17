# ==============================================================================
# 03_relative_trend_distributions.R
# Distribution of relative trends (masked vs unmasked)
# ==============================================================================
suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(here)
})

# ---- config ------------------------------------------------------------------
VAR <- "LAI"
MASK <- "CCI"
TAU <- "tau_0.1"
METRIC <- "yearmean"

EPS_MU <- as.numeric(if (VAR == "LAI") {
  "0.05"
} else {
  "0.02"
})

LAB_UNM <- "Unmasked (post-abiotic)"
LAB_MSK <- sprintf("Natural-only (%s-based)", MASK)
OUT_DIR <- here("analysis", "results", "relative_trend_distribution", TAU)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

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

# ---- paths -------------------------------------------------------------------
AREA025 <- here("src", "area_0p25_validdomain_km2.nc")

REL_UNM <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_%s_trend_relative_peryear_0p25.nc", VAR, METRIC)
)
REL_MSK <- here(
  "output", TAU, "eval",
  sprintf("trend_%s_%s", VAR, MASK),
  sprintf("%s_%s_trend_relative_peryear_0p25.nc", VAR, METRIC)
)

MU_UNM <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_%s_0p25.nc", VAR, METRIC)
)
MU_MSK <- here(
  "output", TAU, "eval",
  sprintf("trend_%s_%s", VAR, MASK),
  sprintf("%s_%s_0p25.nc", VAR, METRIC)
)

# ---- helper ------------------------------------------------------------------
sym_lim <- function(x, q = 0.995, fallback = 2, min_n = 50L) {
  x <- x[is.finite(x)]
  if (length(x) < min_n) {
    return(c(-fallback, fallback))
  }
  lim <- as.numeric(stats::quantile(abs(x), probs = q, na.rm = TRUE))
  if (!is.finite(lim) || lim <= 0) {
    lim <- fallback
  }
  c(-lim, lim)
}

# ---- load (single layer) -----------------------------------------------------
A <- rast(AREA025)[[1]]
r_unm <- rast(REL_UNM)[[1]]
r_msk <- rast(REL_MSK)[[1]]
mu_unm <- rast(MU_UNM)[[1]]
mu_msk <- rast(MU_MSK)[[1]]

# ---- extract values ----------------------------------------------------------
vA <- values(A, mat = FALSE)
vRu <- values(r_unm, mat = FALSE)
vRm <- values(r_msk, mat = FALSE)
vMu <- values(mu_unm, mat = FALSE)
vMm <- values(mu_msk, mat = FALSE)

valid_dom <- is.finite(vA) & (vA > 0)

# filter: require finite trend AND finite mean level above EPS
ok_unm <- valid_dom & is.finite(vRu) & is.finite(vMu) & (vMu > EPS_MU)
ok_msk <- valid_dom & is.finite(vRm) & is.finite(vMm) & (vMm > EPS_MU)

# Use percent per year consistently
df <- bind_rows(
  tibble(Domain = LAB_UNM, r_pct = 100 * vRu[ok_unm], area = vA[ok_unm]),
  tibble(Domain = LAB_MSK, r_pct = 100 * vRm[ok_msk], area = vA[ok_msk])
) |>
  group_by(Domain) |>
  mutate(w = area / sum(area, na.rm = TRUE)) |>
  ungroup()

# ---- summary table -----------------------------------------------------------
dom_area_km2 <- sum(vA[valid_dom], na.rm = TRUE)

sum_tbl <- df |>
  group_by(Domain) |>
  summarise(
    area_valid_km2 = sum(area, na.rm = TRUE),
    area_valid_frac_of_validdom = area_valid_km2 / dom_area_km2,
    reltrend_mean_area_pct_per_year = sum(
      r_pct * area,
      na.rm = TRUE
    ) / sum(area, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(across(where(is.double), ~ round(.x, 4)))

out_csv <- file.path(
  OUT_DIR,
  sprintf(
    "reltrend_distribution_%s_%s_%s_%s_summary.csv", VAR, MASK, TAU, METRIC
  )
)
write_csv(sum_tbl, out_csv)

# ---- figure ------------------------------------------------------------------
lims <- sym_lim(df$r_pct, q = 0.995, fallback = 2)

fill_vals <- setNames(c("grey35", "#2C7FB8"), c(LAB_UNM, LAB_MSK))

p <- ggplot(df, aes(x = r_pct, weight = w, fill = Domain)) +
  geom_histogram(
    bins = 80,
    position = "identity", alpha = 0.55, colour = "white", linewidth = 0.2
  ) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey60") +
  coord_cartesian(xlim = lims) +
  scale_fill_manual(values = fill_vals) +
  labs(
    title = sprintf(
      "%s %s: area-weighted distribution of relative trends", VAR, METRIC
    ),
    subtitle = sprintf(
      "Mask: %s; τ: %s; mean-level filter: %s > %.3f", MASK, TAU, VAR, EPS_MU
    ),
    x = "Relative trend (% yr-1)",
    y = "Area fraction per bin (within domain)",
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank()) +
  theme_pub()

out_png <- file.path(
  OUT_DIR,
  sprintf("reltrend_distribution_%s_%s_%s_%s.png", VAR, MASK, TAU, METRIC)
)
ggsave(out_png, p, width = 11, height = 4.2, dpi = 300)
