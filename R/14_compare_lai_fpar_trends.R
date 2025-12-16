## =============================================================================
## 14_compare_lai_fpar_trends.R — Compare LAI vs FPAR trend structures (0.25°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(scico)
  library(here)
  library(tidyr)
  library(rnaturalearth)
  library(patchwork)
  library(readr)
  library(scales)
})
source(here("R", "utils.R"))
source(here("R", "viz.R"))
ROOT   <- here::here()
DIR025 <- file.path(ROOT, "analysis/unmasked/0p25")
OUTDIR <- file.path(ROOT, "analysis/lai_vs_fpar")

OUTDIR_Q <- file.path(OUTDIR, "quadrants")
dir.create(OUTDIR_Q, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

f_lai  <- file.path(DIR025, "LAI_georef_yearmean_trend_slope_peryear_0p25.nc")
f_fpar <- file.path(DIR025, "FPAR_georef_yearmean_trend_slope_peryear_0p25.nc")
f_lai_base  <- file.path(DIR025, "LAI_georef_yearmean_0p25.nc")
f_fpar_base <- file.path(DIR025, "FPAR_georef_yearmean_0p25.nc")

BASE_START <- 1982
BASE_END   <- 2000
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

# ----------------------------------------------------------------------
# 1. Load rasters efficiently
# ----------------------------------------------------------------------
# Load baseline time series
r_lai_base  <- rast(f_lai_base)
r_fpar_base <- rast(f_fpar_base)

years_base <- BASE_START:(BASE_START + nlyr(r_lai_base) - 1)
idx_base   <- which(years_base >= BASE_START & years_base <= BASE_END)

lai_ref  <- mean(r_lai_base[[idx_base]],  na.rm = TRUE)
fpar_ref <- mean(r_fpar_base[[idx_base]], na.rm = TRUE)

lai_ref  <- ifel(abs(lai_ref)  < 1e-8, NA_real_, lai_ref)
fpar_ref <- ifel(abs(fpar_ref) < 1e-8, NA_real_, fpar_ref)

tol <- 1e-6
ok_ref <- is.finite(lai_ref) & is.finite(fpar_ref) & (abs(lai_ref) > tol) & (abs(fpar_ref) > tol)

r_lai  <- ifel(ok_ref, 100 * rast(f_lai)  / lai_ref,  NA_real_)
r_fpar <- ifel(ok_ref, 100 * rast(f_fpar) / fpar_ref, NA_real_)

r_stack <- c(r_lai, r_fpar)
names(r_stack) <- c("lai_trend", "fpar_trend")

# Extract only finite pixels
df <- as.data.frame(r_stack, xy = TRUE, na.rm = FALSE) |>
  dplyr::rename(lon = x, lat = y) |>
  dplyr::select(lon, lat, lai_trend, fpar_trend) |>
  dplyr::filter(is.finite(lai_trend), is.finite(fpar_trend))
if (nrow(df) == 0)
  stop("No overlapping valid pixels in LAI–fPAR trends.")

# ----------------------------------------------------------------------
# 2. Statistics
# ----------------------------------------------------------------------
lai  <- df$lai_trend
fpar <- df$fpar_trend

pearson_r   <- cor(lai, fpar, method = "pearson")
spearman_r  <- cor(lai, fpar, method = "spearman")
kendall_tau <- cor(lai, fpar, method = "kendall")

lmfit    <- lm(fpar ~ lai)
slope    <- coef(lmfit)[2]
r2       <- summary(lmfit)$r.squared
rmse     <- sqrt(mean((fpar - predict(lmfit))^2))
N        <- length(lai)
sign_agree <- mean(sign(lai) == sign(fpar)) * 100

stats_label <- sprintf(
  "r = %.2f | ρ = %.2f | τ = %.2f\nR² = %.2f | slope = %.2f | RMSE = %.3f\nSign agreement: %.1f%% | N = %d",
  pearson_r,
  spearman_r,
  kendall_tau,
  r2,
  slope,
  rmse,
  sign_agree,
  N
)

# ----------------------------------------------------------------------
# 3. Scatterplot
# ----------------------------------------------------------------------
p_scatter <- ggplot(df, aes(lai_trend, fpar_trend)) +
  geom_point(alpha = 0.15,
             size = 0.4,
             colour = "#000000") +
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "#0072B2") +
  labs(
    x = "LAI trend (% yr⁻¹)",
    y = "fAPAR trend (% yr⁻¹)",
    title = "Pixel-wise correspondence of LAI and fAPAR trends",
    subtitle = stats_label
  ) +
  theme_pub()

ggsave(
  file.path(OUTDIR, "lai_fpar_trend_scatter.png"),
  p_scatter,
  width = 6,
  height = 5.5,
  dpi = 330
)

# ----------------------------------------------------------------------
# 4. Zonal means
# ----------------------------------------------------------------------
df <- df |> mutate(lat_band = floor(lat))

df_zonal <- df |>
  group_by(lat_band) |>
  summarise(
    lai_mean  = weighted.mean(lai_trend,  w = cos(lat * pi / 180), na.rm = TRUE),
    fpar_mean = weighted.mean(fpar_trend, w = cos(lat * pi / 180), na.rm = TRUE),
    .groups = "drop"
  )

p_zonal <- ggplot(df_zonal, aes(lat_band)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_line(aes(y = lai_mean),
            color = "#0072B2",
            linewidth = 0.8) +
  geom_line(aes(y = fpar_mean),
            color = "#D55E00",
            linewidth = 0.8) +
  scale_x_continuous(labels = lab_deg) +
  labs(
    x = "Latitude",
    y = "Relative trend (% yr⁻¹)",
    title = "Zonal comparison of LAI and fAPAR trends",
    subtitle = "1° zonal means"
  ) +
  theme_pub()

ggsave(
  file.path(OUTDIR, "lai_fpar_trend_zonal.png"),
  p_zonal,
  width = 6,
  height = 5,
  dpi = 330
)

# ----------------------------------------------------------------------
# 5. Trend maps
# ----------------------------------------------------------------------


df_lai_map  <- df |> select(lon, lat, slope = lai_trend)
df_fpar_map <- df |> select(lon, lat, slope = fpar_trend)
p_map_lai  <- plot_map_slope(df_lai_map, "LAI")
p_map_fpar <- plot_map_slope(df_fpar_map, "fAPAR")
ggsave(
  file.path(OUTDIR, "map_lai_trend_peryear.png"),
  p_map_lai,
  width = 7.2,
  height = 3.8,
  dpi = 330
)
ggsave(
  file.path(OUTDIR, "map_fpar_trend_peryear.png"),
  p_map_fpar,
  width = 7.2,
  height = 3.8,
  dpi = 330
)
p_side <- (p_map_lai + p_map_fpar) +
  plot_layout(ncol = 2, guides = "collect")
ggsave(
  filename = file.path(OUTDIR, "map_lai_fpar_trends_side_by_side.png"),
  plot = p_side,
  width = 14,
  height = 4,
  dpi = 330
)

cat("Finished: LAI–FPAR trend comparison.\n")

## =============================================================================
## Additional analysis: LAI–fAPAR quadrants & area fractions
## =============================================================================


sign_eps <- function(x, eps = 0.01){
  ifelse(x > eps, 1L, ifelse(x < -eps, -1L, 0L))
}

df_quad <- df |>
  mutate(
    lat_band = floor(lat),
    weight   = cos(lat * pi / 180),
    s_lai  = sign_eps(lai_trend),
    s_fpar = sign_eps(fpar_trend),
    quadrant = case_when(
      s_lai ==  1 & s_fpar ==  1 ~ "Q1_both_pos",
      s_lai == -1 & s_fpar == -1 ~ "Q2_both_neg",
      s_lai ==  1 & s_fpar == -1 ~ "Q3_LAIpos_FPARneg",
      s_lai == -1 & s_fpar ==  1 ~ "Q4_LAIneg_FPARpos",
      TRUE                       ~ "Q0_weak_or_zero"
    )
  )

# ----------------------------------------------------------------------
# 1. Global area fractions
# ----------------------------------------------------------------------
quad_global <- df_quad |>
  group_by(quadrant) |>
  summarise(w_total = sum(weight, na.rm = TRUE), .groups = "drop") |>
  mutate(frac_area = w_total / sum(w_total))

write_csv(quad_global,
          file.path(OUTDIR_Q, "quadrant_global_fractions.csv"))

# ----------------------------------------------------------------------
# 2. Zonal area fractions
# ----------------------------------------------------------------------
quad_zonal <- df_quad |>
  group_by(lat_band, quadrant) |>
  summarise(w_total = sum(weight, na.rm = TRUE), .groups = "drop") |>
  group_by(lat_band) |>
  mutate(frac_area = w_total / sum(w_total)) |>
  ungroup()

write_csv(quad_zonal,
          file.path(OUTDIR_Q, "quadrant_zonal_fractions.csv"))

# ----------------------------------------------------------------------
# 3. Quadrant map
# ----------------------------------------------------------------------

quad_cols <- c(
  Q1_both_pos         = "#0072B2",
  Q2_both_neg         = "#D55E00",
  Q3_LAIpos_FPARneg   = "#009E73",
  Q4_LAIneg_FPARpos   = "#CC79A7",
  Q0_weak_or_zero     = "grey70"
)

p_quad_map <- ggplot(df_quad, aes(lon, lat)) +
  geom_raster(aes(fill = quadrant)) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  coord_sf(expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 60), labels = lab_deg) +
  scale_y_continuous(breaks = seq(-90, 90, 30), labels = lab_deg) +
  scale_fill_manual(values = quad_cols, name   = "LAI–fAPAR trend quadrant") +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "LAI–fAPAR trend quadrants (per year)",
    subtitle = "Q1: both greening; Q2: both browning; Q3/Q4: mixed-sign; Q0: weak/near-zero"
  ) +
  theme_pub()

ggsave(
  file.path(OUTDIR_Q, "quadrant_map.png"),
  p_quad_map,
  width = 7.5,
  height = 4,
  dpi = 330
)

cat("Saved LAI–fAPAR quadrant analysis to:\n  ", OUTDIR_Q, "\n")

