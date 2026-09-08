## =============================================================================
# 03_cci_mask_0p05.R — Build conservative CCI/C3S “used-land” mask (0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))
cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.9)

tmpl <- rast(cfg$grids$grid_005$ref_raster)
out_dir <- cfg$paths$masks_cci_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Mask definition
alpha_cci <- cfg$esa_cci$used_land$threshold
k_cci <- as.integer(cfg$esa_cci$used_land$persistence_years)

cci_years <- cfg$project$years$cci_start:cfg$project$years$cci_end
band_name <- "frac_fused"

frac_dir <- cfg$paths$cci_out_dir
years <- cci_years
fpaths <- file.path(frac_dir, sprintf("ESACCI_frac_%d_0p05.tif", years))
missing_files <- fpaths[!file.exists(fpaths)]
if (length(missing_files)) {
  stop("Missing expected CCI fraction files:\n", paste(missing_files, collapse = "\n"))
}

#  load stack
cci_stack <- rast(lapply(fpaths, function(f) {
  # Load the raster and ensure the band name is correct
  r <- rast(f)
  if (!band_name %in% names(r)) {
    if (nlayers(r) == 2L) {
      names(r) <- c("frac_fused", "frac_grass")
    } else {
      stop(
        sprintf(
          "Band '%s' not found in %s; available bands: %s",
          band_name,
          basename(f),
          paste(names(r), collapse = ", ")
        )
      )
    }
  }
  r[[band_name]]
}))
time(cci_stack) <- years

cci_stack <- align_to_template(cci_stack, tmpl, method = "bilinear")

nl <- nlyr(cci_stack)
k_eff <- min(k_cci, nl)

message(
  sprintf(
    "CCI stack: band=%s, nlayers=%d, years=[%d..%d], alpha=%.3f, k=%d",
    band_name,
    nl,
    min(years),
    max(years),
    alpha_cci,
    k_eff
  )
)

#  majority mask
used_year <- cci_stack >= alpha_cci
mask_log <- app(used_year, sum, na.rm = TRUE) >= k_eff
mask_byte <- ifel(mask_log, 1L, 0L)

y1 <- min(years)
y2 <- max(years)
alpha_tok <- gsub("\\.", "p", sprintf("%.2f", alpha_cci))

out_mask_cci <- file.path(
  out_dir,
  sprintf(
    "mask_used_%s_alpha%s_k%d_%d-%d_0p05.tif",
    band_name,
    alpha_tok,
    k_eff,
    y1,
    y2
  )
)

if (!file.exists(out_mask_cci)) {
  writeRaster(
    mask_byte,
    out_mask_cci,
    overwrite = TRUE,
    wopt = wopt_byte(FALSE)
  )
}

gc()
