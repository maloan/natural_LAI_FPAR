# 11_eval_outputs.R
library(tidyverse)
library(scales)
# --- config & refs -------------------------------------------------------------
# Root dir
ROOT <- normalizePath(path.expand(
  Sys.getenv("SNU_LAI_FPAR_ROOT", "~/GitHub/natural_LAI_FPAR")
),
winslash = "/",
mustWork = FALSE)
# Other source files
source(file.path(ROOT, "R", "00_utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))
source(file.path(ROOT, "R", "viz.R"))
cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)
OUTDIR <- file.path(ROOT, "output", "eval")
dir.create(OUTDIR, TRUE, showWarnings= FALSE)

# 1) Global masked fraction (bar)
read_csv(file.path(OUTDIR, "global_masked_fraction.csv")) %>%
  mutate(mask = factor(mask, levels=c("CCI","GLC","LUH"))) %>%
  ggplot(aes(mask, p_drop, fill = mask)) +
  geom_col(width=0.6) +
  scale_y_continuous(labels = percent_format(accuracy=0.1), expand=expansion(mult=c(0,0.05))) +
  scale_fill_manual(values=c(CCI="#4daf4a", GLC="#377eb8", LUH="#e41a1c")) +
  labs(title="Area fraction removed by each mask",
       x=NULL, y="Area removed (%)") +
  theme_minimal(base_size = 12) + theme(legend.position="none")
ggsave(file.path(OUTDIR, "fig_masked_fraction.png"), width=5.5, height=3.6, dpi=200)

# 2) Pairwise Jaccard + overlap table heatmap
pair <- read_csv(file.path(OUTDIR, "pairwise_comparison_area_weighted.csv"))

# Jaccard bars
pair %>%
  mutate(pair = fct_relevel(pair, "CCI_vs_GLC","CCI_vs_LUH","GLC_vs_LUH")) %>%
  ggplot(aes(pair, jacc, fill=pair)) +
  geom_col(width=0.6) +
  scale_y_continuous(labels=percent_format(accuracy=0.1), limits=c(0,1)) +
  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb")) +
  labs(title="Agreement on dropped area (Jaccard index)", x=NULL, y="Jaccard (drop=1)") +
  theme_minimal(base_size=12) + theme(legend.position="none")
ggsave(file.path(OUTDIR, "fig_pairwise_jaccard.png"), width=6, height=3.6, dpi=200)

# Optional: normalized contingency heatmap for one pair (e.g., CCI vs GLC)
mk_heat <- function(df, which_pair, outfile){
  d <- df %>% filter(pair==which_pair) %>% select(w00,w01,w10,w11)
  mat <- matrix(c(d$w00, d$w01, d$w10, d$w11), nrow=2, byrow=TRUE,
                dimnames=list(CCI=c("keep","drop"), Other=c("keep","drop")))
  matp <- mat / sum(mat)
  tibble(
    cci = rep(c("keep","drop"), each=2),
    other = rep(c("keep","drop"), times=2),
    p = as.vector(matp)
  ) %>%
    ggplot(aes(other, cci, fill=p)) +
    geom_tile() +
    geom_text(aes(label=scales::percent(p, accuracy=0.1)), color="white", size=4) +
    scale_fill_viridis_c(labels=percent_format(accuracy=1)) +
    labs(title=paste0(which_pair,": area-weighted contingency"),
         x="Other mask", y="CCI") +
    theme_minimal(base_size=12) + theme(legend.position="right")
  ggsave(outfile, width=5, height=4, dpi=200)
}
mk_heat(pair, "CCI_vs_GLC", file.path(OUTDIR,"fig_pair_heat_CCI_vs_GLC.png"))


# 3) Sensitivity: CCI τ/k
cci <- read_csv(file.path(OUTDIR, "sensitivity_cci_tau_k.csv"))
cci %>%
  mutate(k = factor(k)) %>%
  ggplot(aes(tau, p_drop, color=k)) +
  geom_line() + geom_point(size=1.6) +
  scale_y_continuous(labels=percent_format(accuracy=0.1)) +
  scale_color_brewer(palette="Set1", name="k (years)") +
  labs(title="CCI used-mask sensitivity", x=expression(tau), y="Area removed (%)") +
  theme_minimal(base_size=12)
ggsave(file.path(OUTDIR, "fig_sens_cci_tau_k.png"), width=6.2, height=3.8, dpi=200)

# 4) Sensitivity: GLC N
glc <- read_csv(file.path(OUTDIR, "sensitivity_glc_N.csv"))
glc %>%
  ggplot(aes(N, p_drop)) +
  geom_line() + geom_point(size=1.8) +
  scale_y_continuous(labels=percent_format(accuracy=0.1)) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(title="GLC used≥N sensitivity", x="N (years)", y="Area removed (%)") +
  theme_minimal(base_size=12)
ggsave(file.path(OUTDIR, "fig_sens_glc_N.png"), width=6.2, height=3.8, dpi=200)

# 5) Sensitivity: LUH α, Gmin, Pmin (heatmap at fixed Gmin or Pmin)
luh <- read_csv(file.path(OUTDIR, "sensitivity_luh_alpha_gmin_pmin.csv"))
# Example: heatmap vs alpha & Pmin at Gmin=0.05
luh %>% filter(abs(gmin - 0.05) < 1e-6) %>%
  ggplot(aes(alpha, pmin, fill=p_drop)) +
  geom_tile(color="white") +
  scale_fill_viridis_c(labels=percent_format(accuracy=0.1), name="Area removed") +
  scale_x_continuous(breaks=unique(sort(luh$alpha))) +
  scale_y_continuous(breaks=unique(sort(luh$pmin))) +
  labs(title="LUH pasture-overlap sensitivity (Gmin=0.05)",
       x=expression(alpha~~"(pasture / grass)"), y=expression(P[min])) +
  theme_minimal(base_size=12)
ggsave(file.path(OUTDIR, "fig_sens_luh_heatmap.png"), width=6.2, height=4.6, dpi=200)


