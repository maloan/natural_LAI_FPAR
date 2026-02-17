## =============================================================================
# 04_cci_mask_0p05.R — Build conservative CCI/C3S “used-land” mask (0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "viz.R"))
cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.9)

tmpl <- rast(cfg$grids$grid_005$ref_raster)
out_dir <- cfg$paths$masks_cci_dir
ql_dir <- file.path(out_dir, "quicklooks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

# --- parameters (env-only defaults) -------------------------------------------
tau_cci <- as.numeric(Sys.getenv("TAU_CCI", "0.1"))
k_cci <- as.integer(Sys.getenv("K_CCI", "3"))

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))

cci_years <- cfg$project$years$cci_start:cfg$project$years$cci_end

# --- files --------------------------------------------------------------------
frac_dir <- cfg$paths$cci_out_dir
fpaths <- list.files(
  frac_dir,
  pattern = "ESACCI_frac_\\d{4}_0p05\\.tif$", full.names = TRUE
)
stopifnot(length(fpaths) > 0)
years <- as.integer(sub("^.*_(\\d{4})_0p05\\.tif$", "\\1", basename(fpaths)))
ord <- order(years)
fpaths <- fpaths[ord]
years <- years[ord]

keep <- years %in% cci_years
fpaths <- fpaths[keep]
years <- years[keep]
stopifnot(length(fpaths) > 0)
# --- validate band ------------------------------------------------------------
r0 <- rast(fpaths[1])

# --- load stack ---------------------------------------------------------------
cci_stack <- rast(lapply(fpaths, function(f) rast(f)[["frac_fused"]]))
time(cci_stack) <- years

if (!compareGeom(cci_stack, tmpl, stopOnError = FALSE)) {
  cci_stack <- resample(cci_stack, tmpl, method = "bilinear")
}

nl <- nlyr(cci_stack)
k_eff <- min(k_cci, nl)

message(sprintf(
  "CCI stack: band=%s, nlayers=%d, years=[%d..%d], tau=%.3f, k=%d",
  "frac_fused", nl, min(years), max(years), tau_cci, k_eff
))

# --- persistence mask ---------------------------------------------------------
used_year <- cci_stack >= tau_cci
mask_log <- app(used_year, sum, na.rm = TRUE) >= k_eff
mask_byte <- ifel(mask_log, 1L, 0L)

y1 <- min(years)
y2 <- max(years)
tau_tok <- gsub("\\.", "p", sprintf("%.2f", tau_cci))
tag <- sprintf("%s_tau%s_k%d_%d-%d", "frac_fused", tau_tok, k_eff, y1, y2)

out_mask_cci <- file.path(out_dir, sprintf(
  "mask_used_%s_tau%s_k%d_%d-%d_0p05.tif",
  "frac_fused", tau_tok, k_eff, y1, y2
))

if (!(SKIP_EXISTING && file.exists(out_mask_cci))) {
  writeRaster(mask_byte, out_mask_cci,
    overwrite = TRUE,
    wopt = wopt_byte(opts$SPEED_OVER_SIZE)
  )
}

# --- quicklooks ---------------------------------------------------------------
ql_probe <- file.path(
  ql_dir, "global",
  sprintf("quicklook_mask_global_%s.png", tag)
)
if (REMAKE_QL || !file.exists(ql_probe)) {
  quicklook_mask_all_aois(
    mask = mask_byte,
    title = sprintf(
      "CCI mask (%s, \u03C4=%.2f, k=%d)",
      "frac_fused", tau_cci, k_eff
    ),
    tag = tag,
    cfg = cfg,
    ql_root = ql_dir,
    down = 4L,
    include_global = TRUE,
    drop_global_key = FALSE
  )
}

gc()
cat(sprintf("Written: %s\nQuicklooks in: %s\n", out_mask_cci, ql_dir))
