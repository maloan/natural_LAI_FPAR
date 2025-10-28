## =============================================================================
# 04_cci_mask_0p05.R — Build conservative CCI/C3S “used-land” mask (0.05°)
#
# Purpose
#   Create a binary “drop” mask (1=drop, 0=keep) from annual ESA-CCI/C3S
#   fractional covers by identifying pixels persistently above a fraction
#   threshold τ across ≥k years within the chosen window.
#
# Inputs
#   - Fraction rasters: ESACCI_frac_{YYYY}_0p05.tif (from 02_cci_frac_0p05.R)
#   - Reference grid: cfg$grids$grid_005$ref_raster
#   - Config thresholds: tau_cci, k_cci, cci_band, cci_years
#
# Outputs
#   - mask_used_{band}_tau{τ}_k{k}_{Y1–Y2}_0p05.tif (Byte; 1=drop, 0=keep, NA=255)
#   - Quicklook PNGs for global and AOIs
#
# Environment variables
#   CCI_BAND       (string, default "frac_fused")
#   TAU_CCI        (numeric, default 0.1)
#   K_CCI          (integer, default 3)
#   SKIP_EXISTING  (logical, default TRUE)
#   REMAKE_QL      (logical, default FALSE)
#
# Dependencies
#   Packages: terra, stringr, glue
#   Sourced helpers: utils.R, io.R, geom.R, viz.R, options.R
#
# Processing overview
#   1) Stack yearly fractional cover rasters (CCI/C3S, 1992–2020).
#   2) Threshold each year by τ to mark used pixels.
#   3) Require persistence ≥k years → conservative “used-land” mask.
#   4) Write Byte GeoTIFF (1=drop, 0=keep) and generate AOI quicklooks.
## =============================================================================



# load packages
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
source(file.path(ROOT, "R", "utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))
source(file.path(ROOT, "R", "viz.R"))
source(file.path(ROOT, "R", "options.R"))
cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)

tmpl    <- rast(cfg$grids$grid_005$ref_raster)
out_dir <- path.expand(cfg$paths$masks_cci_dir)
ql_dir  <- file.path(out_dir, "quicklooks")

dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)
# --- helpers ------------------------------------------------------------------
get2 <- function(lst, nm, default)
  if (!is.null(lst[[nm]])) {
    lst[[nm]]
  } else {
    default
  }

# --- thresholds / band selection ----------------------------------------------
mask_cfg  <- get2(cfg, "mask", list())
cci_band <- Sys.getenv("CCI_BAND", unset = get2(mask_cfg, "cci_band", "frac_fused"))
tau_cci   <- as.numeric(Sys.getenv("TAU_CCI", get2(mask_cfg, "tau_cci", 0.1)))
k_cci     <- as.integer(Sys.getenv("K_CCI", get2(mask_cfg, "k_cci", 3L)))
SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
REMAKE_QL     <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))

cci_years <- get2(mask_cfg, "cci_years", 1992:2020)

# Discover files + build stack for the requested band
frac_dir <- file.path(cfg$paths$cci_out_dir)
fpaths <- list.files(frac_dir, "ESACCI_frac_\\d{4}_0p05\\.tif$", full.names = TRUE)
stopifnot(length(fpaths) > 0)
# Validate band; prefer frac_fused if present
r0 <- rast(fpaths[1])
if (!(cci_band %in% names(r0))) {
  if ("frac_fused" %in% names(r0)) {
    message("Band '", cci_band, "' not found; using 'frac_fused'.")
    cci_band <- "frac_fused"
  } else {
    stop("Neither requested band nor 'frac_fused' found in: ", fpaths[1])
  }
}
# extract 4-digit year from each filename
years <- as.integer(stringr::str_extract(basename(fpaths), "\\d{4}"))
ord   <- order(years)
fpaths <- fpaths[ord] # reindex fpaths/years: sort files by year
years <- years[ord] # reindex fpaths/years: sort files by year
keep   <- years %in% cci_years # keep only years requested in config
fpaths <- fpaths[keep] # Subset fpaths and years to keep
years <- years[keep] # Subset fpaths and years to keep
stopifnot(length(fpaths) > 0) # error if nothing remains.
y1_cfg <- cfg$project$years$cci_start
y2_cfg <- cfg$project$years$cci_end
stopifnot(min(years) >= y1_cfg, max(years) <= y2_cfg)
# Opens each yearly file, extracts the layer named cci_band (e.g., frac_fused),
# and stacks them into a multi-layer SpatRaster (one layer per year).
cci_stack <- rast(lapply(fpaths, function(f)
  rast(f)[[cci_band]]))
# Attach the years vector as the time index for the stack
time(cci_stack) <- years
if (!compareGeom(cci_stack, tmpl, stopOnError = FALSE)) {
  # If grid/CRS/resolution don’t match the template tmpl, resample the stack to
  # that grid using bilinear interpolation
  cci_stack <- resample(cci_stack, tmpl, method = "bilinear")
}

message(
  sprintf(
    "CCI stack: band=%s, nlayers=%d, years=[%d..%d], tau=%.3f, k=%d",
    cci_band,
    nlyr(cci_stack),
    min(years),
    max(years),
    tau_cci,
    k_cci
  )
)

# Build mask (used ≥ k years where band ≥ tau)   -> 1=drop, 0=keep For each
# year/layer, mark pixels that meet the fraction threshold 𝜏 τ. Result: logical
# stack (TRUE/FALSE per year).
used_year <- cci_stack >= tau_cci
nl        <- nlyr(used_year) # Number of years (layers).
k_eff     <- min(k_cci, nl) # Effective persistence requirement: cap𝑘k at the number of available years.
mask_log  <- app(used_year, sum, na.rm = TRUE) >= k_eff # Count, per pixel, in how many years it was TRUE; set TRUE where the count ≥𝑘 eff k eff
mask_byte <- ifel(mask_log, 1L, 0L) # Convert logical mask to byte 0/1 (by your convention: 1=drop, 0=keep).


# Outputs (embed band, tau, k, year window)
y1 <- min(years)
y2 <- max(years)
tau_tok <- gsub("\\.", "p", sprintf("%.2f", tau_cci))
tag <- sprintf("%s_tau%s_k%d_%d-%d", cci_band, tau_tok, k_eff, y1, y2)

out_mask_cci <- file.path(
  out_dir,
  sprintf(
    "mask_used_%s_tau%s_k%d_%d-%d_0p05.tif",
    cci_band,
    tau_tok,
    k_eff,
    y1,
    y2
  )
)

# Write (respect SKIP_EXISTING)
if (!(SKIP_EXISTING && file.exists(out_mask_cci))) {
  writeRaster(
    mask_byte,
    out_mask_cci,
    overwrite = TRUE,
    gdal = gdal_wopt("LOG1S")$gdal,
    NAflag = 255
  )
}

# Quicklooks: Global + all AOIs (folders created automatically)
ql_probe <- file.path(ql_dir,
                      "global",
                      sprintf("quicklook_mask_global_%s.png", tag))

if (REMAKE_QL || !file.exists(ql_probe)) {
  quicklook_mask_all_aois(
    mask    = mask_byte,
    title   = sprintf("CCI mask (%s, \u03C4=%.2f, k=%d)", cci_band, tau_cci, k_eff),
    tag     = tag,
    cfg     = cfg,
    ql_root = ql_dir,
    down    = 4L,
    include_global  = TRUE,
    drop_global_key = FALSE  # set TRUE if cfg$aois includes a 'global' entry
  )
}
gc()
cat(glue("Written: {out_mask_cci}, Quicklooks in: {ql_dir}"))
