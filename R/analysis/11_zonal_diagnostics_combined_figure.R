# ==============================================================================
# 11_zonal_diagnostics_combined_figure.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(here)
})

source(here("R", "helpers", "plotting.R"))

lat_label_fn <- get("lat_labels", mode = "function")
theme_pub_fn <- get("theme_pub", mode = "function")

alpha <- "alpha_0.1"
var <- "LAI"

f_abs <- here(
  "analysis",
  "results",
  "tables",
  "zonal",
  sprintf("zonal_absolute_trends_all_masks_%s.csv", alpha)
)
f_amp <- here(
  "analysis",
  "results",
  "tables",
  "zonal",
  sprintf(
    "zonal_yearamp_timeMean_%s_all_masks_%s.csv",
    toupper(var),
    alpha
  )
)

if (!file.exists(f_abs)) {
  stop("Missing absolute-trends table: ", f_abs)
}
if (!file.exists(f_amp)) {
  stop("Missing seasonal-amplitude table: ", f_amp)
}

df_abs <- read_csv(f_abs, show_col_types = FALSE)
df_amp <- read_csv(f_amp, show_col_types = FALSE)

scenario_levels <- c("Unmasked", "CCI alpha=0.05", "CCI alpha=0.1", "CCI alpha=0.2", "GLC")
df_abs <- df_abs |>
  mutate(scenario = factor(.data$scenario, levels = scenario_levels))
df_amp <- df_amp |>
  mutate(scenario = factor(.data$scenario, levels = scenario_levels))

df_abs_mean <- df_abs |> filter(metric == "Annual mean")
df_abs_max <- df_abs |> filter(metric == "Annual maximum")
abs_ylim <- range(
  c(
    df_abs_mean$abstrend_m2m2yr,
    df_abs_max$abstrend_m2m2yr
  ),
  na.rm = TRUE
)
abs_pad <- 0.05 * diff(abs_ylim)
abs_ylim <- abs_ylim + c(-abs_pad, abs_pad)

p1 <- plot_zonal_diagnostics(
  df_abs_mean,
  "abstrend_m2m2yr",
  "Annual Mean Absolute Trend",
  expression("Absolute trend (% yr"^
    {
      -1
    } * ")"),
  y_limits = abs_ylim
)

p2 <- plot_zonal_diagnostics(
  df_abs_max,
  "abstrend_m2m2yr",
  "Annual Maximum Absolute Trend",
  expression("Absolute trend (% yr"^
    {
      -1
    } * ")"),
  y_limits = abs_ylim
)
p3 <- plot_zonal_diagnostics(
  df_amp,
  "mean_yearamp",
  "Seasonal Amplitude",
  expression("Mean seasonal amplitude LAI (" * m^2 ~ m^
    {
      -2
    } * ")")
)

# Combine panels
fig <- (p1 / p2 / p3) +
  plot_layout(guides = "collect", heights = c(1, 1, 1.1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(plot.caption = element_text(
    size = 9,
    hjust = 0,
    margin = margin(t = 10)
  ))
outdir <- here("analysis", "results", "figures", "summaries")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_stub <- sprintf(
  "zonal_diagnostics_mean_max_amplitude_%s_all_masks_%s",
  toupper(var),
  alpha
)
ggsave(
  file.path(outdir, paste0(out_stub, ".png")),
  fig,
  width = 10.5,
  height = 11.5,
  dpi = 400
)
ggsave(
  file.path(outdir, paste0(out_stub, ".pdf")),
  fig,
  width = 10.5,
  height = 11.5,
  device = cairo_pdf
)
