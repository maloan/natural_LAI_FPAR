## =============================================================================
# 04_cci_mask_0p05.R — Build conservative CCI/C3S “used-land” mask (0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(stringr)
  library(glue)
  library(here)
})

# --- config & refs -------------------------------------------------------------
ROOT <- here()

source(here("R", "utils.R"))
source(here("R", "io.R"))
source(here("R", "geom.R"))
source(here("R", "viz.R"))
source(here("R", "options.R"))

cfg  <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

tmpl    <- rast(cfg$grids$grid_005$ref_raster)
out_dir <- cfg$paths$masks_cci_dir
ql_dir  <- file.path(out_dir, "quicklooks")

dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)

# --- helpers ------------------------------------------------------------------
get2 <- function(lst, nm, default)
  if (!is.null(lst[[nm]])) lst[[nm]] else default

# --- thresholds / band selection ----------------------------------------------
mask_cfg   <- get2(cfg, "mask", list())
cci_band   <- Sys.getenv("CCI_BAND", unset = get2(mask_cfg, "cci_band", "frac_fused"))
tau_cci    <- as.numeric(Sys.getenv("TAU_CCI", get2(mask_cfg, "tau_cci", 0.1)))
k_cci      <- as.integer(Sys.getenv("K_CCI",  get2(mask_cfg, "k_cci", 3L)))
SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
REMAKE_QL     <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))

cci_years <- get2(mask_cfg, "cci_years", 1992:2020)

# Discover files for selected band
frac_dir <- cfg$paths$cci_out_dir
fpaths <- list.files(frac_dir, "ESACCI_frac_\\d{4}_0p05\\.tif$", full.names = TRUE)
stopifnot(length(fpaths) > 0)

# Validate band
r0 <- rast(fpaths[1])
if (!(cci_band %in% names(r0))) {
  if ("frac_fused" %in% names(r0)) {
    message("Band '", cci_band, "' not found; using 'frac_fused'.")
    cci_band <- "frac_fused"
  } else {
    stop("Neither requested band nor 'frac_fused' found in: ", fpaths[1])
  }
}

years <- as.integer(stringr::str_extract(basename(fpaths), "\\d{4}"))
ord   <- order(years)
fpaths <- fpaths[ord]
years  <- years[ord]

keep <- years %in% cci_years
fpaths <- fpaths[keep]
years  <- years[keep]
stopifnot(length(fpaths) > 0)

y1_cfg <- cfg$project$years$cci_start
y2_cfg <- cfg$project$years$cci_end
stopifnot(min(years) >= y1_cfg, max(years) <= y2_cfg)

cci_stack <- rast(lapply(fpaths, function(f) rast(f)[[cci_band]]))
time(cci_stack) <- years

if (!compareGeom(cci_stack, tmpl, stopOnError = FALSE)) {
  cci_stack <- resample(cci_stack, tmpl, method = "bilinear")
}

message(
  sprintf(
    "CCI stack: band=%s, nlayers=%d, years=[%d..%d], tau=%.3f, k=%d",
    cci_band, nlyr(cci_stack),
    min(years), max(years),
    tau_cci, k_cci
  )
)

# Build mask (persistence ≥ k)
used_year <- cci_stack >= tau_cci
nl        <- nlyr(used_year)
k_eff     <- min(k_cci, nl)

mask_log  <- app(used_year, sum, na.rm = TRUE) >= k_eff
mask_byte <- ifel(mask_log, 1L, 0L)

# --- output filenames ----------------------------------------------------------
y1 <- min(years)
y2 <- max(years)

tau_tok <- gsub("\\.", "p", sprintf("%.2f", tau_cci))
tag <- sprintf("%s_tau%s_k%d_%d-%d", cci_band, tau_tok, k_eff, y1, y2)

out_mask_cci <- file.path(
  out_dir,
  sprintf(
    "mask_used_%s_tau%s_k%d_%d-%d_0p05.tif",
    cci_band, tau_tok, k_eff, y1, y2
  )
)

# Write mask
if (!(SKIP_EXISTING && file.exists(out_mask_cci))) {
  writeRaster(
    mask_byte,
    out_mask_cci,
    overwrite = TRUE,
    gdal = gdal_wopt("LOG1S")$gdal,
    NAflag = 255
  )
}

# --- quicklooks ----------------------------------------------------------------
ql_probe <- file.path(
  ql_dir,
  "global",
  sprintf("quicklook_mask_global_%s.png", tag)
)

if (REMAKE_QL || !file.exists(ql_probe)) {
  quicklook_mask_all_aois(
    mask    = mask_byte,
    title   = sprintf("CCI mask (%s, \u03C4=%.2f, k=%d)", cci_band, tau_cci, k_eff),
    tag     = tag,
    cfg     = cfg,
    ql_root = ql_dir,
    down    = 4L,
    include_global  = TRUE,
    drop_global_key = FALSE
  )
}

gc()
cat(glue("Written: {out_mask_cci}, Quicklooks in: {ql_dir}"))
