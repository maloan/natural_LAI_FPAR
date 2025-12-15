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

ROOT   <- here::here()
DIR025 <- file.path(ROOT, "analysis/unmasked/0p25")
OUTDIR <- file.path(ROOT, "analysis/lai_vs_fpar")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

f_lai  <- file.path(DIR025, "LAI_georef_yearmean_trend_slope_peryear_0p25.nc")
f_fpar <- file.path(DIR025, "FPAR_georef_yearmean_trend_slope_peryear_0p25.nc")

# ----------------------------------------------------------------------
# Theme
# ----------------------------------------------------------------------
theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey87", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 13, face = "bold"),
      plot.subtitle    = element_text(size = 10),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 9)
    )
}

lab_deg <- scales::label_number(suffix = "°")
coast   <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

# ----------------------------------------------------------------------
# 1. Load rasters efficiently
# ----------------------------------------------------------------------
r_lai  <- rast(f_lai)
r_fpar <- rast(f_fpar)

# Stack without building a huge df
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
intercept <- coef(lmfit)[1]
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
    x = "LAI trend (per year)",
    y = "fAPAR trend (per year)",
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
df_lai_map  <- df |> select(lon, lat, slope = lai_trend)
df_fpar_map <- df |> select(lon, lat, slope = fpar_trend)
df <- df |> mutate(lat_band = floor(lat))

df_zonal <- df |>
  group_by(lat_band) |>
  summarise(
    lai_mean  = mean(lai_trend, na.rm = TRUE),
    fpar_mean = mean(fpar_trend, na.rm = TRUE),
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
    y = "Trend (per year)",
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
plot_trend_map_fast <- function(df, var_name, SD_K = 5) {
  sdev  <- sd(df$slope, na.rm = TRUE)
  clamp <- 2 * sdev
  df$slope_clamped <- pmax(pmin(df$slope, clamp), -clamp)
  ggplot(df, aes(lon, lat)) +
    geom_raster(aes(fill = slope_clamped)) +
    geom_sf(
      data = coast,
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.2
    ) +
    coord_sf(expand = FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 60), labels = lab_deg) +
    scale_y_continuous(breaks = seq(-90, 90, 30), labels = lab_deg) +
    scale_fill_scico(
      palette = "bam",
      name    = paste0(var_name, " trend (per year)"),
      limits  = c(-clamp, clamp),
      oob     = scales::squish
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title    = paste0(var_name, " trend (per year)")
    ) +
    theme_pub()
}

df_lai_map  <- df |> select(lon, lat, slope = lai_trend)
df_fpar_map <- df |> select(lon, lat, slope = fpar_trend)
p_map_lai  <- plot_trend_map_fast(df_lai_map, "LAI")
p_map_fpar <- plot_trend_map_fast(df_fpar_map, "fAPAR")
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

OUTDIR_Q <- file.path(OUTDIR, "quadrants")
dir.create(OUTDIR_Q, recursive = TRUE, showWarnings = FALSE)

EPS_Q <- 0  # threshold for classifying sign; change if needed (e.g. 1e-4)

sign_eps <- function(x, eps = 0) {
  ifelse(x > eps, 1L, ifelse(x < -eps, -1L, 0L))
}

df_quad <- df |>
  mutate(
    s_lai  = sign_eps(lai_trend, EPS_Q),
    s_fpar = sign_eps(fpar_trend, EPS_Q),
    quadrant = case_when(
      s_lai ==  1 & s_fpar ==  1 ~ "Q1_both_pos",
      s_lai == -1 & s_fpar == -1 ~ "Q2_both_neg",
      s_lai ==  1 & s_fpar == -1 ~ "Q3_LAIpos_FPARneg",
      s_lai == -1 & s_fpar ==  1 ~ "Q4_LAIneg_FPARpos",
      TRUE                       ~ "Q0_weak_or_zero"
    ),
    weight   = cos(lat * pi / 180),
    lat_band = floor(lat)
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
coast   <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
lab_deg <- scales::label_number(suffix = "°")

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

# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------


ROOT   <- here::here()
INDIR  <- file.path(ROOT, "analysis", "lai_vs_fpar", "quadrants")
OUTDIR <- file.path(INDIR, "plots")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(file.path(INDIR, "quadrant_zonal_fractions.csv"))

theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey87", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold")
    )
}

lab_deg <- scales::label_number(suffix = "°")

quad_cols <- c(
  Q1_both_pos         = "#1b9e77",
  Q2_both_neg         = "#d95f02",
  Q3_LAIpos_FPARneg   = "#7570b3",
  Q4_LAIneg_FPARpos   = "#e7298a",
  Q0_weak_or_zero     = "grey60"
)

p <- ggplot(df, aes(lat_band, frac_area, colour = quadrant)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = quad_cols, name = NULL) +
  scale_x_continuous(labels = lab_deg) +
  labs(x = "Latitude (°)", y = "Area fraction", title = "Zonal fractions of LAI–fAPAR trend quadrants") +
  theme_pub() +
  theme(legend.position = "bottom")

ggsave(
  file.path(OUTDIR, "lai_fpar_quadrant_zonal.png"),
  p,
  width = 7,
  height = 5,
  dpi = 330
)

cat("Saved LAI–fAPAR quadrant zonal plot to:\n  ", OUTDIR, "\n")

