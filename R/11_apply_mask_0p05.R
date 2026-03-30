## =============================================================================
# 11_apply_mask_0p05.R — Apply drop masks to monthly LAI/FPAR at 0.05°
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "viz.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)

var <- toupper(Sys.getenv("var", "FPAR"))
mask <- toupper(Sys.getenv("mask", "CCI"))

skip_existing <- as_bool(Sys.getenv("skip_existing"), default = TRUE)
overwrite <- as_bool(Sys.getenv("overwrite"), default = FALSE)
src <- toupper(Sys.getenv("grass_source", "CCI"))
g_min <- as.numeric(Sys.getenv("g_min", "0.1"))
p_min <- as.numeric(Sys.getenv("p_min", "0.1"))
alpha <- as.numeric(Sys.getenv("alpha", "0.50"))

ref005 <- rast(cfg$grids$grid_005$ref_raster)

in_dir <- if (var == "LAI") {
  cfg$paths$georef_lai_0p05_dir
} else {
  cfg$paths$georef_fpar_0p05_dir
}
out_dir <- switch(mask,
  CCI = if (var == "LAI") {
    cfg$paths$masked_lai_cci_0p05_dir
  } else {
    cfg$paths$masked_fpar_cci_0p05_dir
  },
  GLC = if (var == "LAI") {
    cfg$paths$masked_lai_glc_0p05_dir
  } else {
    cfg$paths$masked_fpar_glc_0p05_dir
  }
)

stopifnot(is.character(out_dir), length(out_dir) == 1, nzchar(out_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

year_1 <- cfg$project$years$cci_start
year_2 <- cfg$project$years$cci_end

# Format parameter values for filename matching
tok <- function(x) gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))

mask_path <- if (mask == "CCI") {
  find_one(
    cfg$paths$masks_cci_dir,
    sprintf("mask_used_.*_%d-%d_0p05\\.tif$", year_1, year_2)
  )
} else {
  num <- as.integer(Sys.getenv("USED_N_YEARS", "3"))
  find_one(
    cfg$paths$masks_glc_dir,
    sprintf("mask_used_ge%d_%d-%d_0p05\\.tif$", num, year_1, year_2)
  )
}

mask <- rast(mask_path)
if (!compareGeom(mask, ref005, stopOnError = FALSE)) {
  mask <- resample(mask, ref005, method = "near")
}

combine_or <- function(a, b) {
  if (!compareGeom(b, a, stopOnError = FALSE)) {
    b <- resample(b, a, method = "near")
  }
  app(c(a, b), fun = function(v) as.integer(any(v >= 1, na.rm = TRUE)))
}

abi_dir <- file.path(cfg$paths$masks_root_dir, "mask_abiotic")
abi_path <- list.files(abi_dir, "mask_abiotic_CCI_2007_.*_0p05\\.tif$",
  full.names = TRUE
)
if (length(abi_path)) {
  abi <- rast(abi_path[order(file.info(abi_path)$mtime, decreasing = TRUE)][1])
  mask <- combine_or(mask, abi)
}

rx <- sprintf(
  "mask_luh_overlap_%s_Gmin%s_Pmin%s_alpha%s_%d-%d_0p05_rep\\.tif$",
  src, tok(g_min), tok(p_min), tok(alpha), year_1, year_2
)
luh_path <- find_one(
  file.path(cfg$paths$masks_root_dir, "mask_luh_overlap"), rx
)
luh <- rast(luh_path)
mask <- combine_or(mask, luh)
vals_ok <- try(all(values(mask) %in% c(0, 1, NA)), silent = TRUE)

mask_combined_path <- file.path(out_dir, "combined_mask_0p05.tif")
if (!file.exists(mask_combined_path) || overwrite) {
  writeRaster(mask, mask_combined_path,
    overwrite = TRUE,
    wopt = wopt_byte(opts$speed_over_size, na = 255L)
  )
}

files <- sort(list.files(
  in_dir,
  pattern = paste0("^", var, "_\\d{6}_0p05\\.tif$"), full.names = TRUE
))
stopifnot(length(files) > 0L)

wopt <- wopt_f32(opts$speed_over_size)

for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_%s_0p05_masked.tif", var, ym))

  do_write <- overwrite || !file.exists(out) || !skip_existing

  r <- rast(f)

  if (do_write) {
    r_masked <- mask(r, mask, maskvalues = 1, updatevalue = NA)
    writeRaster(r_masked, out, overwrite = TRUE, wopt = wopt)
  }
}

gc()
message("Masked ", var, " written to: ", out_dir)
