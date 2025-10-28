## =============================================================================
# 07_apply_mask_0p05.R — Apply drop masks to monthly LAI/FPAR at 0.05°

suppressPackageStartupMessages({
  library(terra)
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

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "FALSE"))
OVERWRITE     <- as.logical(Sys.getenv("OVERWRITE", "TRUE"))
REMAKE_QL     <- as.logical(Sys.getenv("REMAKE_QL", "TRUE"))
APPLY_LUH_OVERLAP <- as.logical(Sys.getenv("APPLY_LUH_OVERLAP","TRUE"))
APPLY_ABIOTIC_STATIC <- as.logical(Sys.getenv("APPLY_ABIOTIC_STATIC", "TRUE"))
Y0   <- as.integer(Sys.getenv("LUH_AVG_START", cfg$project$years$cci_start))
YL   <- as.integer(Sys.getenv("LUH_AVG_END",   cfg$project$years$cci_end))
SRC  <- toupper(Sys.getenv("GRASS_SOURCE", "CCI"))   # CCI|GLC_FRAC|GLC_TEMP
GMIN <- as.numeric(Sys.getenv("G_MIN",  "0.05"))
PMIN <- as.numeric(Sys.getenv("P_MIN",  "0.20"))

find_one <- function(dir, pattern) {
  cand <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(cand))
    stop("No mask matching pattern '", pattern, "' in: ", dir, call. = FALSE)
  cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]  # newest
}

# --- choose variable + mask source --------------------------------------------
VAR  <- toupper(Sys.getenv("VAR", "FPAR"))     # LAI | FPAR
MASK <- toupper(Sys.getenv("MASK", "CCI"))    # CCI | GLC | LUH_RATIO
USED_N_YEARS <- as.integer(Sys.getenv("USED_N_YEARS", "3"))
stopifnot(VAR %in% c("LAI", "FPAR"), MASK %in% c("CCI", "GLC", "LUH_RATIO"))

# --- io, clamps, patterns ------------------------------------------------------
ref005 <- rast(cfg$grids$grid_005$ref_raster)

if (VAR == "LAI") {
  in_dir   <- cfg$paths$georef_lai_0p05_dir
  clamp_lo <- cfg$variables$lai$clamp$min
  clamp_hi <- cfg$variables$lai$clamp$max
  patt     <- "^LAI_\\d{6}_0p05\\.tif$"
  ql_title <- "LAI"
  zlim     <- c(clamp_lo, clamp_hi)
  out_tmpl_cci <- "LAI_SNUv1_cu_tau0p15_k3_{ym}_0p05_masked.tif"
  out_tmpl_glc <- glue("LAI_SNUv1_glc_ge{USED_N_YEARS}_{ '{ym}' }_0p05_masked.tif")
  out_tmpl_luh <- "LAI_SNUv1_luh_ratio_{ym}_0p05_masked.tif"
  out_dir <- switch(
    MASK,
    "CCI"       = cfg$paths$masked_lai_cci_005_dir,
    "GLC"       = cfg$paths$masked_lai_glc_005_dir,
    "LUH_RATIO" = cfg$paths$masked_lai_luh_005_dir
  )
} else {
  in_dir   <- cfg$paths$georef_fpar_0p05_dir
  clamp_lo <- cfg$variables$fpar$clamp$min
  clamp_hi <- cfg$variables$fpar$clamp$max
  patt     <- "^FPAR_\\d{6}_0p05\\.tif$"
  ql_title <- "FPAR"
  zlim     <- c(clamp_lo, clamp_hi)
  out_tmpl_cci <- "FPAR_SNUv1_cu_tau0p15_k3_{ym}_0p05_masked.tif"
  out_tmpl_glc <- glue("FPAR_SNUv1_glc_ge{USED_N_YEARS}_{ '{ym}' }_0p05_masked.tif")
  out_tmpl_luh <- "FPAR_SNUv1_luh_ratio_{ym}_0p05_masked.tif"
  out_dir <- switch(
    MASK,
    "CCI"       = cfg$paths$masked_fpar_cci_005_dir,
    "GLC"       = cfg$paths$masked_fpar_glc_005_dir,
    "LUH_RATIO" = cfg$paths$masked_fpar_luh_005_dir
  )
}
out_tmpl <- switch(MASK,
                   "CCI" = out_tmpl_cci,
                   "GLC" = out_tmpl_glc,
                   "LUH_RATIO" = out_tmpl_luh)

# Ensure output dirs exist
dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "quicklooks"), TRUE, showWarnings = FALSE)

# --- pick mask file (semantics: 1 = drop, 0 = keep) ---------------------------
Y1 <- cfg$project$years$cci_start
Y2 <- cfg$project$years$cci_end

if (MASK == "CCI") {
  band_env <- Sys.getenv("CCI_BAND", "")  # e.g., frac_fused
  tau_env  <- Sys.getenv("TAU_CCI", "")
  k_env    <- Sys.getenv("K_CCI", "")
  tau_tok  <- if (nzchar(tau_env))
{    gsub("\\.", "p", sprintf("%.2f", as.numeric(tau_env)))}
  else
    {".*"}

  if (nzchar(tau_env) && nzchar(k_env)) {
    rx <- if (nzchar(band_env)) {
      sprintf("mask_used_%s_tau%s_k%s_%d-%d_0p05\\.tif$",
              band_env,
              tau_tok,
              k_env,
              Y1,
              Y2)
    } else {
      sprintf("mask_used_.*_tau%s_k%s_%d-%d_0p05\\.tif$",
              tau_tok,
              k_env,
              Y1,
              Y2)
    }
  } else {
    rx <- if (nzchar(band_env)) {
      sprintf("mask_used_%s_tau.*_k\\d+_%d-%d_0p05\\.tif$",
              band_env,
              Y1,
              Y2)
    } else {
      sprintf("mask_used_.*_tau.*_k\\d+_%d-%d_0p05\\.tif$", Y1, Y2)
    }
  }
  mask_path <- find_one(cfg$paths$masks_cci_dir, rx)
  # sanity: mask year window matches config (embedded in filename)
  rxwin <- "(\\d{4})-(\\d{4})"
  m <- regmatches(mask_path, regexec(rxwin, mask_path))[[1]]
  if (length(m) == 3) {
    stopifnot(as.integer(m[2]) == Y1, as.integer(m[3]) == Y2)
  }

} else if (MASK == "GLC") {
  N <- as.integer(Sys.getenv("USED_N_YEARS", USED_N_YEARS))
  rx <- sprintf("mask_used_ge%d_%d-%d_0p05\\.tif$", N, Y1, Y2)
  mask_path <- find_one(cfg$paths$masks_glc_dir, rx)

}

message("Using mask: ", mask_path)

# Align mask to template (categorical → nearest)
mask <- rast(mask_path)
if (!compareGeom(mask, ref005, stopOnError = FALSE)) {
  mask <- resample(mask, ref005, method = "near")
}
if (!all(values(mask) %in% c(0,1,NA))) stop("Unexpected mask values — expected 0/1 only.")

mask_combined_path <- file.path(out_dir, "combined_mask_0p05.tif")
if (!file.exists(mask_combined_path) || OVERWRITE) {
  writeRaster(mask, mask_combined_path, overwrite=TRUE,
              gdal = gdal_wopt("LOG1S")$gdal, NAflag=255)
}
# --- discover inputs -----------------------------------------------------------
files <- list.files(in_dir, patt, full.names = TRUE)
if (!length(files))
  stop("No ", VAR, " files found in: ", in_dir)

# --- quicklook helper ----------------------------------------------------------
quicklook_after_full <- function(ra, ym, title = ql_title, down = 4L, zlim = zlim) {
  qd <- file.path(out_dir, "quicklooks")
  dir.create(qd, TRUE, showWarnings = FALSE)
  rr <- if (down > 1L) aggregate(ra, down, mean, na.rm = TRUE) else ra
  png(file.path(qd, sprintf("quicklook_%s_masked_full_%s.png", title, ym)),
      width = 1400, height = 700, res = 120)
  op <- par(oma = c(0,0,2,0), mar = c(3,3,3,6))
  terra::plot(rr,
              main  = sprintf("%s %s (masked)", title, ym),
              col   = pal_green(64), colNA = col_na, zlim = zlim,
              axes  = TRUE, legend = TRUE
  )
  .add_overlays(rr)
  mtext("Longitude (°E)", side = 1, line = 2)
  mtext("Latitude (°N)",  side = 2, line = 2)
  par(op); dev.off()
}


quicklook <- function(rb,
                      ra,
                      ym,
                      title = ql_title,
                      down = 4L,
                      zlim = zlim) {
  qd <- file.path(out_dir, "quicklooks")
  dir.create(qd, TRUE, showWarnings = FALSE)
  if (down > 1L) {
    rb <- aggregate(rb, down, mean, na.rm = TRUE)
    ra <- aggregate(ra, down, mean, na.rm = TRUE)
  }
  png(
    file.path(qd, sprintf(
      "quicklook_%s_masked_%s.png", title, ym
    )),
    width = 1400,
    height = 700,
    res = 120
  )
  op <- par(
    mfrow = c(1, 2),
    oma = c(0, 0, 2.2, 0),
    mar = c(3, 3, 3, 6)
  )
  terra::plot(
    rb,
    main = sprintf("%s %s (before)", title, ym),
    col = pal_green(64),
    colNA = col_na,
    zlim = zlim,
    axes = TRUE,
    legend = TRUE
  )
  .add_overlays(rb)
  mtext("Longitude (°E)", side = 1, line = 2)
  mtext("Latitude (°N)", side = 2, line = 2)
  terra::plot(
    ra,
    main = sprintf("%s %s (masked)", title, ym),
    col = pal_green(64),
    colNA = col_na,
    zlim = zlim,
    axes = TRUE,
    legend = TRUE
  )
  .add_overlays(ra)
  mtext("Longitude (°E)", side = 1, line = 2)
  mtext("Latitude (°N)", side = 2, line = 2)
  mtext(
    "Masked 0.05° quicklook",
    side = 3,
    outer = TRUE,
    cex = 1.05
  )
  par(op)
  dev.off()
}

# --- optional static abiotic overlay (union) ----------------------------------

if (APPLY_ABIOTIC_STATIC) {
  abi_dir  <- file.path(cfg$paths$masks_root_dir, "mask_abiotic")
  abi_path <- list.files(abi_dir, "mask_static_abiotic_CCI_.*_0p05.tif", full.names = TRUE)
   if (length(abi_path)) {
    abi_path <- abi_path[order(file.info(abi_path)$mtime, decreasing = TRUE)][1]
    abi_mask <- rast(abi_path)
    if (!compareGeom(abi_mask, ref005, stopOnError = FALSE)) {
      abi_mask <- resample(abi_mask, ref005, method = "near")
    }
    mask <- app(c(mask, abi_mask),
                fun = function(v) as.integer(any(v >= 1, na.rm = TRUE)),
                cores = 1)

    message("Applied static abiotic overlay: ", basename(abi_path))
  } else {
    message("APPLY_ABIOTIC_STATIC=TRUE, but no abiotic mask found in ",
            abi_dir)
  }
}

# --- optional LUH overlay (pasture-overlap) -----------------------------------
if (APPLY_LUH_OVERLAP) {
  tok <- function(x) gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
  rx <- glue("mask_luh_overlap_CCI_Gmin{tok(GMIN)}_Pmin{tok(PMIN)}_alpha{tok(ALPHA)}_{Y1}-{Y2}_0p05_rep.tif")
  luh_dir  <- file.path(cfg$paths$masks_root_dir, "mask_luh_overlap")
  luh_path <- find_one(luh_dir, rx)
  luh <- rast(luh_path)
  if (!compareGeom(luh, ref005, stopOnError = FALSE)) {
    luh <- resample(luh, ref005, method = "near")
  }
  # union: 1=drop if either mask says drop
  mask <- app(c(mask, luh),
              fun = function(v) as.integer(any(v >= 1, na.rm = TRUE)),
              cores = 1)
  message("Applied LUH pasture-overlap overlay: ", basename(luh_path))
}

# --- clamp + mask loop ---------------------------------------------------------
for (f in files) {
  ym <- str_extract(basename(f), "\\d{6}")
  out <- file.path(out_dir, glue(out_tmpl))

  if (!(SKIP_EXISTING && file.exists(out)) || OVERWRITE) {
    r <- rast(f)
    if (!compareGeom(r, ref005, stopOnError = FALSE)) {
      stop("Geometry mismatch for ",
           basename(f),
           " vs ref005; regrid upstream.")
    }
    r <- clamp(r, lower = clamp_lo, upper = clamp_hi)
    r_masked <- mask(r, mask, maskvalues = 1, updatevalue = NA)  # 1=drop, 0=keep
    # Assert clamp still respected
    mm <- minmax(r_masked)
    stopifnot(mm[1, 1] >= clamp_lo - 1e-6 || is.na(mm[1, 1]))
    stopifnot(mm[2, 1] <= clamp_hi + 1e-6 || is.na(mm[2, 1]))

    writeRaster(
      r_masked,
      out,
      overwrite = TRUE,
      gdal   = gdal_wopt("FLT4S")$gdal,
      NAflag = -9999
    )
  } else {
    # still load for QL if needed
    r <- rast(f)
    r_masked <- rast(out)
  }


  # Jan/Jul quicklooks
  if (substr(ym, 5, 6) %in% c("01", "07")) {
    ql_png <- file.path(out_dir,"quicklooks",
                        sprintf("quicklook_%s_masked_%s.png", ql_title, ym))
    if (REMAKE_QL || !(SKIP_EXISTING && file.exists(ql_png))) {
      quicklook(r, r_masked, ym, title = ql_title, zlim = zlim)
    }
    # after-only, full-frame
    ql_full <- file.path(out_dir,"quicklooks",
                         sprintf("quicklook_%s_masked_full_%s.png", ql_title, ym))
    if (REMAKE_QL || !(SKIP_EXISTING && file.exists(ql_full))) {
      quicklook_after_full(r_masked, ym, title = ql_title, zlim = zlim)
    }
  }

}

message("Masked ", VAR, " written to: ", out_dir)
