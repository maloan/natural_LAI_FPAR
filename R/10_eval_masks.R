#!/usr/bin/env Rscript
# 10_eval_masks.R — QA, stats, cross-compare, and sensitivity for masks
# Minimal deps: terra, yaml, stringr, glue

suppressPackageStartupMessages({
  library(terra)
  library(yaml)
  library(stringr)
  library(glue)
})
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
stopifnot(!is.null(cfg$paths), !is.null(cfg$grids))

ref005 <- rast(cfg$grids$grid_005$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)  # km^2
stopifnot(!is.null(area005))

# Output folder
OUTDIR <- file.path(ROOT, "output", "eval"); dir.create(OUTDIR, TRUE, showWarnings= FALSE)
QLDIR  <- file.path(OUTDIR, "quicklooks"); dir.create(QLDIR, TRUE, showWarnings= FALSE)

# Helpers
read_mask_safe <- function(p) {
  r <- rast(p); if (!compareGeom(r, ref005, stopOnError = FALSE)) {
    r <- resample(r, ref005, method = "near")  # masks are categorical/binary
  }
  r
}
aw_mean <- function(x, w) {  # area-weighted mean (ignores NA in x and w)
  num <- global(x * w, "sum", na.rm = TRUE)[1,1]
  den <- global((!is.na(x)) * w, "sum", na.rm = TRUE)[1,1]
  if (is.na(den) || den == 0) NA_real_ else as.numeric(num / den)
}
p_drop <- function(mask01, w = area005) {
  # mask01 semantics: 1=drop, 0=keep; returns area-weighted fraction dropped
  aw_mean(ifel(mask01 >= 1, 1, 0), w)
}
pairwise_stats <- function(a, b, w = area005) {
  # area-weighted contingency + Jaccard (on "drop")
  A <- ifel(a >= 1, 1, 0); B <- ifel(b >= 1, 1, 0)
  wA <- global(A * w, "sum", na.rm = TRUE)[1,1]
  wB <- global(B * w, "sum", na.rm = TRUE)[1,1]
  wAND <- global((A & B) * w, "sum", na.rm = TRUE)[1,1]
  wOR  <- global(((A | B)) * w, "sum", na.rm = TRUE)[1,1]
  w00  <- global(((A==0) & (B==0)) * w, "sum", na.rm = TRUE)[1,1]
  w01  <- global(((A==0) & (B==1)) * w, "sum", na.rm = TRUE)[1,1]
  w10  <- global(((A==1) & (B==0)) * w, "sum", na.rm = TRUE)[1,1]
  w11  <- wAND
  jacc <- if (isTRUE(wOR>0)) wAND / wOR else NA_real_
  data.frame(wA, wB, w00, w01, w10, w11, wAND, wOR, jacc)
}

# ==== locate latest masks (adjust patterns if needed) ==========================
# CCI mask (0.05°): mask_used_{band}_tau{...}_k{...}_{Y1}-{Y2}_0p05.tif
cci_dir <- cfg$paths$masks_cci_dir
cci_files <- list.files(cci_dir, "mask_used_frac_fused_tau0p10_k3_1992-2020_0p05.tif", full.names = TRUE)
stopifnot(length(cci_files) > 0)
cci_path <- cci_files[order(file.info(cci_files)$mtime, decreasing = TRUE)][1]

# GLC mask (0.05°): mask_used_ge{N}_{Y1}-{Y2}_0p05.tif
glc_dir <- cfg$paths$masks_glc_dir
glc_files <- list.files(glc_dir, "mask_used_ge3_1992-2020_0p05.tif", full.names = TRUE)
stopifnot(length(glc_files) > 0)
glc_path <- glc_files[order(file.info(glc_files)$mtime, decreasing = TRUE)][1]

# LUH overlap mask (0.05° replica); be flexible with naming:
# e.g., mask_luh_overlap_*_0p05.tif OR mask_pasture_ratio_*_0p05.tif
luh_dir <- cfg$paths$masks_luh_dir
luh_files <- list.files(luh_dir, "m_LUH_pasture_1992-2020_0p05_rep.tif", full.names = TRUE, ignore.case = TRUE)
stopifnot(length(luh_files) > 0)
luh_path <- luh_files[order(file.info(luh_files)$mtime, decreasing = TRUE)][1]

message("Using masks:\n  CCI: ", basename(cci_path),
        "\n  GLC: ", basename(glc_path),
        "\n  LUH: ", basename(luh_path))

m_cci <- read_mask_safe(cci_path)
m_glc <- read_mask_safe(glc_path)
m_luh <- read_mask_safe(luh_path)

# ==== 1) Visual QA quicklooks (per-AOI) =======================================
plot_mask <- function(r, title, file, pal=c("#f0f0f0","#d73027")) {
  png(file, 1400, 700, res=120)
  op <- par(mar=c(3,3,3,6))
  terra::plot(r, main=title, col=pal, breaks=c(-0.5,0.5,1.5),
              legend=FALSE, axes=TRUE, box=TRUE)
  legend("bottomleft", fill=pal, legend=c("0 keep","1 drop"), bty="n")
  par(op); dev.off()
}

# Global quicklooks
plot_mask(m_cci, "CCI used-mask (1=drop)", file.path(QLDIR, "mask_cci_global.png"))
plot_mask(m_glc, "GLC used≥N mask (1=drop)", file.path(QLDIR, "mask_glc_global.png"))
plot_mask(m_luh, "LUH overlap mask (1=drop)", file.path(QLDIR, "mask_luh_global.png"))

# AOI quicklooks (subset a few to keep it light)
aois <- cfg$aois %||% cfg$project$aois
pick <- intersect(c("europe","corn_belt","se_asia","australia_se","amazon","global"),
                  names(aois))
for (nm in pick) {
  a <- aois[[nm]]
  e <- ext(a$lon_min, a$lon_max, a$lat_min, a$lat_max)
  for (pair in list(
    list(m=m_cci, tag="cci"),
    list(m=m_glc, tag="glc"),
    list(m=m_luh, tag="luh")
  )) {
    r <- crop(pair$m, e, snap="out")
    plot_mask(r, sprintf("%s (%s)", toupper(pair$tag), nm),
              file.path(QLDIR, sprintf("mask_%s_%s.png", pair$tag, nm)))
  }
}

# ==== 2) Fraction of masked area globally (area-weighted) =====================
stats <- data.frame(
  mask = c("CCI","GLC","LUH"),
  p_drop = c(p_drop(m_cci), p_drop(m_glc), p_drop(m_luh))
)
write.csv(stats, file.path(OUTDIR, "global_masked_fraction.csv"), row.names = FALSE)
message(glue("Global masked fraction (area-weighted):\n",
             "  CCI={sprintf('%.3f', stats$p_drop[stats$mask=='CCI'])}\n",
             "  GLC={sprintf('%.3f', stats$p_drop[stats$mask=='GLC'])}\n",
             "  LUH={sprintf('%.3f', stats$p_drop[stats$mask=='LUH'])}"))

# ==== 3) Cross-comparison: CCI vs GLC vs LUH ==================================
cmp_cci_glc <- pairwise_stats(m_cci, m_glc)
cmp_cci_luh <- pairwise_stats(m_cci, m_luh)
cmp_glc_luh <- pairwise_stats(m_glc, m_luh)

cmp <- rbind(
  data.frame(pair="CCI_vs_GLC", cmp_cci_glc),
  data.frame(pair="CCI_vs_LUH", cmp_cci_luh),
  data.frame(pair="GLC_vs_LUH", cmp_glc_luh)
)
write.csv(cmp, file.path(OUTDIR, "pairwise_comparison_area_weighted.csv"), row.names = FALSE)

message("Jaccard (area-weighted, on drop=1): ",
        glue("CCI~GLC={sprintf('%.3f', cmp$jacc[cmp$pair=='CCI_vs_GLC'])}, "),
        glue("CCI~LUH={sprintf('%.3f', cmp$jacc[cmp$pair=='CCI_vs_LUH'])}, "),
        glue("GLC~LUH={sprintf('%.3f', cmp$jacc[cmp$pair=='GLC_vs_LUH'])}"))

# ==== 4) Sensitivity analyses (read masks already produced) ====================
# 4a) CCI: vary τ and k (scan files present)
cci_all <- list.files(cci_dir, "^mask_used_.*_tau.*_k\\d+_\\d{4}-\\d{4}_0p05\\.tif$", full.names = TRUE)
if (length(cci_all)) {
  parse_tau  <- function(s) as.numeric(sub("p",".", str_match(s, "_tau([0-9p]+)_")[,2]))
  parse_k    <- function(s) as.integer(str_match(s, "_k(\\d+)_")[,2])
  df <- data.frame(path=cci_all,
                   tau=parse_tau(cci_all),
                   k=parse_k(cci_all),
                   stringsAsFactors = FALSE)
  df <- df[order(df$tau, df$k),]
  df$p_drop <- NA_real_
  for (i in seq_len(nrow(df))) {
    r <- read_mask_safe(df$path[i]); df$p_drop[i] <- p_drop(r)
  }
  write.csv(df, file.path(OUTDIR, "sensitivity_cci_tau_k.csv"), row.names = FALSE)
}

# 4b) GLC: vary N (scan files present)
glc_all <- list.files(glc_dir, "^mask_used_ge\\d+_\\d{4}-\\d{4}_0p05\\.tif$", full.names = TRUE)
if (length(glc_all)) {
  parse_N <- function(s) as.integer(str_match(s, "mask_used_ge(\\d+)_")[,2])
  dg <- data.frame(path=glc_all, N=parse_N(glc_all), stringsAsFactors = FALSE)
  dg <- dg[order(dg$N),]
  dg$p_drop <- NA_real_
  for (i in seq_len(nrow(dg))) {
    r <- read_mask_safe(dg$path[i]); dg$p_drop[i] <- p_drop(r)
  }
  write.csv(dg, file.path(OUTDIR, "sensitivity_glc_N.csv"), row.names = FALSE)
}

# 4c) LUH: vary α, G_MIN, P_MIN (best-effort parse from filenames)
# Accepts names like "...alpha0p50...", "...Gmin0p05...", "...Pmin0p20..." etc.
luh_all <- unique(c(
  list.files(luh_dir, ".*0p05.*\\.tif$", full.names = TRUE),
  list.files(file.path(cfg$paths$masks_root_dir), ".*0p05.*\\.tif$",
             full.names = TRUE, recursive = TRUE)
))
if (length(luh_all)) {
  ex <- function(re, s) {
    m <- str_match(s, re)[,2]; suppressWarnings(as.numeric(sub("p",".", m)))
  }
  dl <- data.frame(path=luh_all,
                   alpha = ex("alpha([0-9p]+)", luh_all),
                   gmin  = ex("Gmin([0-9p]+)",  luh_all),
                   pmin  = ex("Pmin([0-9p]+)",  luh_all),
                   stringsAsFactors = FALSE)
  dl <- dl[!is.na(dl$alpha) | !is.na(dl$gmin) | !is.na(dl$pmin), ]
  if (nrow(dl)) {
    dl$p_drop <- NA_real_
    for (i in seq_len(nrow(dl))) {
      r <- read_mask_safe(dl$path[i]); dl$p_drop[i] <- p_drop(r)
    }
    dl <- dl[order(dl$alpha, dl$gmin, dl$pmin, na.last = TRUE),]
    write.csv(dl, file.path(OUTDIR, "sensitivity_luh_alpha_gmin_pmin.csv"), row.names = FALSE)
  }
}

for (VAR in c("FPAR", "LAI")) {
  for (MASK_NAME in c("CCI", "GLC")) {
    masked_dir <- cfg$paths[[paste0("masked_", tolower(VAR), "_", tolower(MASK_NAME), "_005_dir")]]
    masked_files <- list.files(masked_dir, pattern = paste0("^", VAR, "_\\d{6}_0p05_masked\\.tif$"), full.names = TRUE)
    stopifnot(length(masked_files) > 0)

    # Load stack of all months
    ras_stack <- rast(masked_files)
    ras_mean <- mean(ras_stack, na.rm = TRUE)  # temporal mean

    # Select mask
    mask_obj <- switch(MASK_NAME,
                       "CCI" = m_cci,
                       "GLC" = m_glc,
                       "LUH" = m_luh,
                       stop("Invalid mask"))

    # Subset rasters
    masked_vals   <- mask(ras_mean, mask_obj, maskvalues = 0, inverse = TRUE)  # where mask == 1 (drop)
    unmasked_vals <- mask(ras_mean, mask_obj, maskvalues = 1, inverse = TRUE)  # where mask == 0 (keep)

    # Area-weighted means
    mean_masked   <- aw_mean(masked_vals, area005)
    mean_unmasked <- aw_mean(unmasked_vals, area005)

    # Output
    result <- data.frame(
      variable = VAR,
      mask = MASK_NAME,
      mean_masked = mean_masked,
      mean_unmasked = mean_unmasked
    )

    write.csv(result, file.path(OUTDIR, paste0("mean_", tolower(VAR), "_", tolower(MASK_NAME), "_masked_vs_unmasked.csv")), row.names = FALSE)

    message(glue("{VAR}, {MASK_NAME}: mean (unmasked) = {sprintf('%.3f', mean_unmasked)}, mean (masked) = {sprintf('%.3f', mean_masked)}"))


  }
}


# ==== Done =====================================================================
cat(glue("
Wrote:
  - Global quicklooks: {QLDIR}/mask_*_global.png
  - AOI quicklooks   : {QLDIR}/mask_*_<AOI>.png
  - Global masked fractions: {OUTDIR}/global_masked_fraction.csv
  - Pairwise comparison     : {OUTDIR}/pairwise_comparison_area_weighted.csv
  - Sensitivity (CCI τ,k)   : {OUTDIR}/sensitivity_cci_tau_k.csv
  - Sensitivity (GLC N)     : {OUTDIR}/sensitivity_glc_N.csv
  - Sensitivity (LUH α,Gmin,Pmin): {OUTDIR}/sensitivity_luh_alpha_gmin_pmin.csv

Notes:
  • Fractions are area-weighted by cfg$grids$grid_005$area_raster (km^2).
  • Mask semantics assumed 1=drop, 0=keep; categorical resampling uses 'near'.
  • LUH file name parsing is best-effort; ensure tokens like alpha0p50, Gmin0p05, Pmin0p20 appear.
"))
