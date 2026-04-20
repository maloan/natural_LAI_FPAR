## =============================================================================
# 03_cci_mask_0p05.R — Build conservative CCI/C3S “used-land” mask (0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "paths.R"))
source(here("R", "helpers", "files.R"))
source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))
cfg <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.9)

tmpl <- rast(cfg$grids$grid_005$ref_raster)
out_dir <- cfg$paths$masks_cci_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- parameters (env-only defaults) -------------------------------------------
tau_cci <- as.numeric(Sys.getenv("tau_cci", "0.1"))
k_cci <- as.integer(Sys.getenv("k_cci", "3"))
skip_existing <- as_bool(Sys.getenv("skip_existing"), default = TRUE)

cci_years <- cfg$project$years$cci_start:cfg$project$years$cci_end
band_name <- "frac_fused"

# --- files --------------------------------------------------------------------
frac_dir <- cfg$paths$cci_out_dir
fpaths <- list.files(
  frac_dir,
  pattern = "ESACCI_frac_\\d{4}_0p05\\.tif$", full.names = TRUE
)
stopifnot(length(fpaths) > 0)

# Extract years, sort, and filter to target range
years <- as.integer(sub("^.*_(\\d{4})_0p05\\.tif$", "\\1", basename(fpaths)))
ord <- order(years)
fpaths <- fpaths[ord]
years <- years[ord]

# Keep only years within config range
keep <- years %in% cci_years
fpaths <- fpaths[keep]
years <- years[keep]
stopifnot(length(fpaths) > 0)

# --- load stack ---------------------------------------------------------------
cci_stack <- rast(lapply(fpaths, function(f) {
  r <- rast(f)
  if (!band_name %in% names(r)) {
    if (nlayers(r) == 2L) {
      names(r) <- c("frac_fused", "frac_grass")
    } else {
      stop(sprintf(
        "Band '%s' not found in %s; available bands: %s",
        band_name, basename(f), paste(names(r), collapse = ", ")
      ))
    }
  }
  r[[band_name]]
}))
time(cci_stack) <- years

cci_stack <- align_to_template(cci_stack, tmpl, method = "bilinear")

nl <- nlyr(cci_stack)
k_eff <- min(k_cci, nl)

message(sprintf(
  "CCI stack: band=%s, nlayers=%d, years=[%d..%d], tau=%.3f, k=%d",
  band_name, nl, min(years), max(years), tau_cci, k_eff
))

# --- persistence mask ---------------------------------------------------------
used_year <- cci_stack >= tau_cci
mask_log <- app(used_year, sum, na.rm = TRUE) >= k_eff
mask_byte <- ifel(mask_log, 1L, 0L)

y1 <- min(years)
y2 <- max(years)
tau_tok <- gsub("\\.", "p", sprintf("%.2f", tau_cci))

out_mask_cci <- file.path(out_dir, sprintf(
  "mask_used_%s_tau%s_k%d_%d-%d_0p05.tif",
  band_name, tau_tok, k_eff, y1, y2
))

if (!(skip_existing && file.exists(out_mask_cci))) {
  writeRaster(mask_byte, out_mask_cci,
    overwrite = TRUE,
    wopt = wopt_byte(opts$speed_over_size)
  )
}

gc()
