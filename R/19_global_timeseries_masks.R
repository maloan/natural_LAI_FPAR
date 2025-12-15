
## =============================================================================
## 19_global_timeseries_masks.R — Global annual mean LAI/FPAR (masked vs unmasked)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(readr)
})

ROOT    <- here::here()
IN_GEO  <- file.path(ROOT, "analysis", "unmasked", "0p25")
OUTDIR  <- file.path(ROOT, "analysis", "global_timeseries_masks")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

VARS    <- c("LAI", "FPAR")
METRICS <- c("yearmean", "yearmax")
MASKS   <- c("CCI", "GLC")
TAU     <- c("tau_0.05", "tau_0.1", "tau_0.2")

# ----------------------------------------------------------------------
# Helper
# ----------------------------------------------------------------------
theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey87", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 13, face = "bold"),
      plot.subtitle    = element_text(size = 10),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 9),
      legend.position  = "bottom"
    )
}

build_label <- function(domain, mask = NA_character_, tau = NA_character_) {
  if (domain == "unmasked") {
    return("Unmasked (all land)")
  }
  sprintf("%s (τ=%s)", mask, sub("tau_", "", tau))
}

# ----------------------------------------------------------------------
# Main loop over VAR × METRIC
# ----------------------------------------------------------------------
for (var in VARS) {
  for (met in METRICS) {
    message("==== ", var, " / ", met, " ====")

    ts_list <- list()

    # -------------------------------
    # 1) Unmasked global time series
    # -------------------------------
    f_geo <- file.path(IN_GEO, sprintf("%s_georef_%s_0p25.nc", var, met))

    if (file.exists(f_geo)) {
      r_geo  <- rast(f_geo)
      years  <- 1982:(1982 + nlyr(r_geo) - 1L)

      df_geo <- terra::global(r_geo, "mean", na.rm = TRUE) |>
        as_tibble() |>
        dplyr::rename(value = 1) |>
        mutate(
          year   = years,
          domain = build_label("unmasked"),
          mask   = NA_character_,
          tau    = NA_character_
        )

      ts_list[[length(ts_list) + 1]] <- df_geo
    } else {
      warning("Missing unmasked file: ", f_geo)
    }

    # -------------------------------
    # 2) Masked global time series
    # -------------------------------
    for (tau in TAU) {
      for (mask in MASKS) {
        f_mask <- file.path(
          ROOT,
          "output",
          tau,
          "eval",
          sprintf("trend_%s_%s", var, mask),
          sprintf("%s_%s_0p25.nc", var, met)
        )

        if (!file.exists(f_mask)) {
          next
        }

        r_mask <- rast(f_mask)
        years  <- 1982:(1982 + nlyr(r_mask) - 1L)

        df_mask <- terra::global(r_mask, "mean", na.rm = TRUE) |>
          as_tibble() |>
          dplyr::rename(value = 1) |>
          mutate(
            year   = years,
            domain = build_label("masked", mask, tau),
            mask   = mask,
            tau    = tau
          )

        ts_list[[length(ts_list) + 1]] <- df_mask
      }
    }

    if (length(ts_list) == 0) {
      warning("No time series found for ", var, " / ", met)
      next
    }

    df_ts <- bind_rows(ts_list)

    # ------------------------------------------------------------------
    # Plot: all domains overlapped
    # ------------------------------------------------------------------
    p <- ggplot(df_ts, aes(year, value, colour = domain)) +
      geom_line(linewidth = 0.6) +
      labs(
        x = "Year",
        y = sprintf("%s (annual global %s)", var, met),
        title = sprintf("%s %s: global annual mean — masked vs unmasked", var, met),
        subtitle = "Unmasked (all land) vs natural-only masks (CCI, various τ)"
      ) +
      theme_pub()

    outf <- file.path(OUTDIR,
                      sprintf("%s_%s_global_timeseries_all_masks.png", var, met))

    ggsave(outf,
           p,
           width = 7.5,
           height = 4.8,
           dpi = 330)

    # Also save the table for later analysis
    write_csv(df_ts, file.path(
      OUTDIR,
      sprintf("%s_%s_global_timeseries_all_masks.csv", var, met)
    ))

    message("Saved: ", outf)
  }
}

cat("Finished global time series (masked vs unmasked).\nOutput in:\n  ",
    OUTDIR,
    "\n")
