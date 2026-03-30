# ==============================================================================
# plot_global_timeseries.R
# Global area-weighted annual time series (unmasked vs masked)
# + linear trend lines (OLS) with slope annotation
# + shaded 95% CI ribbons around annual means (effective-n approximation)
# + adds a second mask case (GLC)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(tibble)
  library(ggplot2)
  library(here)
  library(readr)
})

# ---- config ------------------------------------------------------------------
tau <- "tau_0.1"
var <- "LAI"
metric <- "yearmean"

masks <- c("GLC", "CCI")
outdir <- here("analysis", "results", "paper_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- files -------------------------------------------------------------------
f_unm <- here(
  "analysis", "unmasked", "0p25",
  sprintf("%s_georef_%s_0p25.nc", var, metric)
)

f_msk <- function(mask) {
  here(
    "output", tau, "eval", sprintf("trend_%s_%s", var, mask),
    sprintf("%s_%s_0p25.nc", var, metric)
  )
}

f_area <- here("src", "area_0p25_validdomain_km2.nc")

stopifnot(file.exists(f_unm), file.exists(f_area))
stopifnot(all(vapply(masks, function(m) file.exists(f_msk(m)), logical(1))))

# ---- read data ---------------------------------------------------------------
r_unm <- rast(f_unm)
r_msk <- lapply(masks, function(m) rast(f_msk(m)))
names(r_msk) <- masks
area <- rast(f_area)[[1]]

# geometry checks
compareGeom(r_unm, area, stopOnError = TRUE)
for (m in masks) compareGeom(r_unm, r_msk[[m]], stopOnError = TRUE)

# ---- years -------------------------------------------------------------------
years <- seq_len(nlyr(r_unm)) + 1981
y_lab <- "Leaf Area Index"

# ---- area-weighted annual mean + 95% CI --------------------------------------
global_wmean_ci_stack <- function(x_stack, area) {
  area_valid <- mask(area, area > 0)

  out <- lapply(seq_len(nlyr(x_stack)), function(i) {
    x <- mask(x_stack[[i]], area_valid)

    w <- mask(area_valid, is.finite(x))
    sw <- global(w, "sum", na.rm = TRUE)[1, 1]
    if (!is.finite(sw) || sw <= 0) {
      return(c(
        mean = NA_real_, se = NA_real_,
        lo = NA_real_, hi = NA_real_, n_eff = NA_real_
      ))
    }

    num <- global(x * w, "sum", na.rm = TRUE)[1, 1]
    mu <- as.numeric(num / sw)

    sw2 <- global(w * w, "sum", na.rm = TRUE)[1, 1]
    n_eff <- as.numeric((sw * sw) / sw2)

    v_num <- global((x - mu)^2 * w, "sum", na.rm = TRUE)[1, 1]
    v_w <- as.numeric(v_num / sw)

    if (!is.finite(v_w) || v_w < 0 || !is.finite(n_eff) || n_eff <= 1) {
      se <- NA_real_
    } else {
      se <- sqrt(v_w / n_eff)
    }

    lo <- mu - 1.96 * se
    hi <- mu + 1.96 * se

    c(mean = mu, se = se, lo = lo, hi = hi, n_eff = n_eff)
  })

  as.data.frame(do.call(rbind, out))
}

# compute all series
res_unm <- global_wmean_ci_stack(r_unm, area)

res_mask <- lapply(masks, function(m) global_wmean_ci_stack(r_msk[[m]], area))
names(res_mask) <- masks

# assemble long-ish data frame (wide columns, simple)
df <- tibble(
  year = years,
  unmasked = res_unm$mean,
  unmasked_lo = res_unm$lo,
  unmasked_hi = res_unm$hi
)

for (m in masks) {
  df[[paste0("masked_", m)]] <- res_mask[[m]]$mean
  df[[paste0("masked_", m, "_lo")]] <- res_mask[[m]]$lo
  df[[paste0("masked_", m, "_hi")]] <- res_mask[[m]]$hi
}

# ---- linear trends (OLS) -----------------------------------------------------
fmt_slope <- function(x) formatC(x, format = "f", digits = 2)
fmt_p <- function(p) {
  if (!is.finite(p)) {
    return("NA")
  }
  if (p < 1e-3) "<0.001" else formatC(p, format = "f", digits = 3)
}

fit_unm <- lm(unmasked ~ year, data = df)
slope_unm <- coef(fit_unm)[["year"]] * 100
p_unm <- summary(fit_unm)$coefficients["year", "Pr(>|t|)"]

ann_lines <- c(
  paste0(
    "OLS slope (Unmasked): ", fmt_slope(slope_unm),
    "% yr\u207B\u00B9 (p ", fmt_p(p_unm), ")"
  )
)

for (m in masks) {
  ycol <- paste0("masked_", m)
  fit <- lm(df[[ycol]] ~ df$year)
  slope <- coef(fit)[["df$year"]] * 100
  pval <- summary(fit)$coefficients["df$year", "Pr(>|t|)"]
  ann_lines <- c(
    ann_lines, paste0(
      "OLS slope (Masked-", m, "): ",
      fmt_slope(slope), "% yr\u207B\u00B9 (p ", fmt_p(pval), ")"
    )
  )
}

ann_text <- paste(ann_lines, collapse = "\n")

# ---- theme ---------------------------------------------
theme_paper <- function() {
  theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(size = 9),
      axis.text        = element_text(size = 8),
      plot.title       = element_text(size = 11, face = "bold"),
      plot.subtitle    = element_text(size = 9),
      legend.title     = element_blank(),
      legend.text      = element_text(size = 8),
      legend.position  = "bottom"
    )
}

# ---- plot --------------------------------------------------------------------
col_unm <- "black"
col_msk <- c(CCI = "#D55E00", GLC = "#0072B2")

p <- ggplot(df, aes(x = year)) +
  # ribbons
  geom_ribbon(aes(ymin = unmasked_lo, ymax = unmasked_hi),
    fill = col_unm, alpha = 0.12
  ) +
  (  # add ribbons for each mask
    lapply(masks, function(m) {
      ggplot2::geom_ribbon(
        aes(
          ymin = .data[[paste0("masked_", m, "_lo")]],
          ymax = .data[[paste0("masked_", m, "_hi")]]
        ),
        fill = col_msk[[m]],
        alpha = 0.10
      )
    })
  ) +
  # mean lines
  geom_line(aes(y = unmasked, colour = "Unmasked"), linewidth = 0.6) +
  ( # add masked lines
    lapply(masks, function(m) {
      ggplot2::geom_line(
        aes(y = .data[[paste0("masked_", m)]], colour = paste0("Masked-", m)),
        linewidth = 0.6
      )
    })
  ) +
  # trend lines
  geom_smooth(aes(y = unmasked),
    method = "lm", se = FALSE,
    colour = col_unm, linewidth = 0.3, linetype = "22"
  ) +
  ( # masked trend lines
    lapply(masks, function(m) {
      ggplot2::geom_smooth(
        aes(y = .data[[paste0("masked_", m)]]),
        method = "lm", se = FALSE,
        colour = col_msk[[m]], linewidth = 0.3, linetype = "22"
      )
    })
  ) +
  scale_colour_manual(values = c(
    "Unmasked" = col_unm,
    "Masked-CCI" = col_msk[["CCI"]],
    "Masked-GLC" = col_msk[["GLC"]]
  )) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
  labs(
    x = "Year",
    y = y_lab,
    title = sprintf("%s global annual mean", var),
    subtitle = sprintf(
      "\u03C4 = %s (masks: %s)", tau, paste(masks, collapse = ", ")
    )
  ) +
  annotate(
    "text",
    x = min(df$year) - 1,
    y = max(c(
      df$unmasked_hi,
      df$masked_CCI_hi, df$masked_GLC_hi
    ), na.rm = TRUE),
    label = ann_text,
    hjust = 0, vjust = 0.9,
    size = 2.5,
    colour = "black"
  ) +
  theme_paper()

print(p)

# ---- save outputs -------------------------------------------------------------
out_png <- file.path(
  outdir,
  sprintf("%s_%s_global_timeseries_unmasked_CCI_GLC_%s.png", var, metric, tau)
)
out_csv <- file.path(
  outdir,
  sprintf("%s_%s_global_timeseries_unmasked_CCI_GLC_%s.csv", var, metric, tau)
)

ggsave(out_png, p, width = 6.5, height = 3.8, dpi = 400)
write_csv(df, out_csv)
