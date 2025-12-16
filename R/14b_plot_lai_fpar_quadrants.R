## =============================================================================
## 14b_plot_lai_fpar_quadrants.R
## Zonal fractions of LAI–fAPAR trend quadrants
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(here)
})

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
ROOT   <- here::here()
INDIR  <- file.path(ROOT, "analysis", "lai_vs_fpar", "quadrants")
OUTDIR <- file.path(INDIR, "plots")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Shared visualization utilities
source(here("R", "viz.R"))

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
df <- read_csv(file.path(INDIR, "quadrant_zonal_fractions.csv"),
               show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Ensure complete latitude coverage (no implicit interpolation)
# -----------------------------------------------------------------------------
df_plot <- df |>
  tidyr::complete(
    lat_band = seq(-90, 89),
    quadrant,
    fill = list(frac_area = 0)
  )

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
p <- ggplot(df_plot, aes(lat_band, frac_area, colour = quadrant)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = quad_palette()) +
  scale_x_continuous(labels = lab_deg) +
  labs(
    x = "Latitude (°)",
    y = "Area fraction",
    title = "Zonal fractions of LAI–fAPAR trend quadrants"
  ) +
  theme_pub() +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(OUTDIR, "lai_fpar_quadrant_zonal.png"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 330
)

cat("Saved LAI–fAPAR quadrant zonal plot to:\n  ", OUTDIR, "\n")
