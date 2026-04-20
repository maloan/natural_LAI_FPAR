## =============================================================================
# 10_apply_mask_0p05.R — Apply drop masks to monthly LAI/FPAR at 0.05°
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "paths.R"))
source(here("R", "helpers", "files.R"))
source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "plotting.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)

var <- toupper(Sys.getenv("var", "LAI"))
mask_kind <- toupper(Sys.getenv("mask", "GLC"))
if (!var %in% c("LAI", "FPAR")) {
  stop_msg("Unsupported var: ", var, ". Use LAI or FPAR")
}
if (!mask_kind %in% c("CCI", "GLC")) {
  stop_msg("Unsupported mask: ", mask_kind, ". Use CCI or GLC")
}

skip_existing <- as_bool(Sys.getenv("skip_existing"), default = TRUE)
overwrite <- as_bool(Sys.getenv("overwrite"), default = FALSE)
g_min <- as.numeric(Sys.getenv("g_min", "0.1"))
p_min <- as.numeric(Sys.getenv("p_min", "0.1"))
alpha <- as.numeric(Sys.getenv("alpha", "0.50"))

ref005 <- rast(cfg$grids$grid_005$ref_raster)

in_dir <- if (var == "LAI") {
  cfg$paths$georef_lai_0p05_dir
} else {
  cfg$paths$georef_fpar_0p05_dir
}
out_dir <- switch(mask_kind,
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

year_1 <- if (mask_kind == "CCI") {
  cfg$project$years$cci_start
} else {
  cfg$project$years$glc_start
}
year_2 <- if (mask_kind == "CCI") {
  cfg$project$years$cci_end
} else {
  cfg$project$years$glc_end
}


mask_path <- if (mask_kind == "CCI") {
  find_one(
    cfg$paths$masks_cci_dir,
    sprintf("^mask_used_.*_%d-%d_0p05\\.tif$", year_1, year_2),
    label = "CCI used mask"
  )
} else {
  find_one(
    cfg$paths$masks_glc_dir,
    sprintf("^mask_used_ge[0-9]+_%d-%d_0p05\\.tif$", year_1, year_2),
    label = "GLC used mask"
  )
}

drop_mask <- rast(mask_path)
drop_mask <- align_to_template(drop_mask, ref005, method = "near")
stopifnot(compareGeom(drop_mask, ref005, stopOnError = FALSE))

combine_or <- function(a, b) {
  b <- align_to_template(b, a, method = "near")
  app(c(a, b), fun = function(v) as.integer(any(v >= 1, na.rm = TRUE)))
}

nonveg_dir <- file.path(cfg$paths$masks_root_dir, "mask_nonvegetated")
nonveg_path <- list.files(nonveg_dir, "mask_nonvegetated_CCI_2007_.*_0p05\\.tif$",
  full.names = TRUE
)
nonveg <- rast(nonveg_path[order(file.info(nonveg_path)$mtime, decreasing = TRUE)][1])
drop_mask <- combine_or(drop_mask, nonveg)


luh_dir <- file.path(cfg$paths$masks_root_dir, "mask_luh_overlap")
luh_rx <- sprintf(
  "mask_luh_overlap_%s_Gmin%s_Pmin%s_alpha%s_%d-%d_0p05_rep\\.tif$",
  mask_kind, tok(g_min), tok(p_min), tok(alpha), year_1, year_2
)
luh_path <- find_one(luh_dir, luh_rx, label = "LUH overlap mask")
luh <- rast(luh_path)
drop_mask <- combine_or(drop_mask, luh)
vals_ok <- try(all(values(drop_mask) %in% c(0, 1, NA)), silent = TRUE)
if (inherits(vals_ok, "try-error") || !isTRUE(vals_ok)) {
  stop_msg("Combined mask has values outside {0,1,NA}")
}

mask_combined_path <- file.path(out_dir, "combined_mask_0p05.tif")
if (!file.exists(mask_combined_path) || overwrite) {
  writeRaster(drop_mask, mask_combined_path,
    overwrite = TRUE,
    wopt = wopt_byte(opts$speed_over_size, na = 255L)
  )
}
drop_mask <- rast(here("output", "tau_  ")) # sanity check
plot(drop_mask, main = sprintf("Combined mask (%s)", mask_kind))
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
    r_masked <- terra::mask(r, drop_mask, maskvalues = 1, updatevalue = NA)
    writeRaster(r_masked, out, overwrite = TRUE, wopt = wopt)
  }
}

gc()
message("Masked ", var, " written to: ", out_dir)
