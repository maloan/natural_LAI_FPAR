# ==============================================================================
# 11_zonal_diagnostics_combined_figure.R
# Figure 4: Combined zonal diagnostics (mean trend, max trend, amplitude)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(here)
})

plot_theme_pub <- function(base_size = 12, include_legend = FALSE, include_strip = FALSE) {
  theme_list <- list(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = base_size + 1, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = base_size - 1),
    axis.text = element_text(size = 9)
  )
  if (include_strip) {
    theme_list$strip.text <- element_text(size = 10, face = "bold")
  }
  if (include_legend) {
    theme_list$legend.position <- "bottom"
    theme_list$legend.box <- "vertical"
    theme_list$legend.text <- element_text(size = 9)
  }
  theme_bw(base_size = base_size) + theme(!!!theme_list)
}

theme_pub <- plot_theme_pub

tau <- "tau_0.1"
var <- "LAI"

f_rel <- here(
  "analysis", "results", "figures", "summaries",
  sprintf("zonal_relative_trends_all_masks_%s.csv", tau)
)
f_amp <- here(
  "analysis", "results", "figures", "summaries",
  sprintf("zonal_yearamp_timeMean_%s_all_masks_%s.csv", toupper(var), tau)
)

if (!file.exists(f_rel)) {
  stop("Missing relative-trends table: ", f_rel)
}
if (!file.exists(f_amp)) {
  stop("Missing seasonal-amplitude table: ", f_amp)
}

df_rel <- read_csv(f_rel, show_col_types = FALSE)
df_amp <- read_csv(f_amp, show_col_types = FALSE)

col_map <- c(
  "Unmasked" = "grey20",
  "CCI tau=0.05" = "#1b9e77",
  "CCI tau=0.1" = "#d95f02",
  "CCI tau=0.2" = "#7570b3",
  "GLC" = "#386cb0"
)

lat_labels <- function(x) {
  ifelse(x == 0, "0°", ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N")))
}

mk_panel <- function(df, ycol, ttl, ylab) {
  ggplot(df, aes(.data$lat_band, .data[[ycol]], colour = .data$scenario)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_line(linewidth = 0.75, na.rm = TRUE) +
    scale_colour_manual(values = col_map, drop = FALSE) +
    scale_x_continuous(limits = c(-90, 90), breaks = seq(-90, 90, by = 30), labels = lat_labels) +
    labs(title = ttl, x = "Latitude", y = ylab, colour = NULL) +
    plot_theme_pub(base_size = 11, include_legend = TRUE)
}

p1 <- mk_panel(
  df_rel |> filter(metric == "Annual mean"),
  "reltrend_pct_per_year",
  "Annual mean relative trend",
  "% yr-1"
)

p2 <- mk_panel(
  df_rel |> filter(metric == "Annual maximum"),
  "reltrend_pct_per_year",
  "Annual maximum relative trend",
  "% yr-1"
)

p3 <- mk_panel(
  df_amp,
  "mean_yearamp",
  "Seasonal amplitude",
  sprintf("%s amplitude", toupper(var))
)

fig <- (p1 / p2 / p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Latitudinal diagnostics across masking scenarios",
    subtitle = "Rows: annual mean trends, annual maximum trends, and seasonal amplitude"
  ) &
  theme(legend.position = "bottom")

outdir <- here("analysis", "results", "figures", "summaries")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_stub <- sprintf(
  "zonal_diagnostics_mean_max_amplitude_%s_all_masks_%s",
  toupper(var),
  tau
)

ggsave(file.path(outdir, paste0(out_stub, ".png")), fig, width = 10.5, height = 11.5, dpi = 320)
ggsave(file.path(outdir, paste0(out_stub, ".pdf")), fig, width = 10.5, height = 11.5)
