# =============================================================================
# 10_apply_mask_0p05.R — Apply drop masks to monthly LAI/FPAR at 0.05°
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "plotting.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()
terraOptions(progress = 1, memfrac = 0.25)

var <- toupper(Sys.getenv("VAR", "LAI"))
mask_kind <- toupper(Sys.getenv("MASK", "GLC"))
if (!var %in% c("LAI", "FPAR")) {
  stop_msg("Unsupported var: ", var, ". Use LAI or FPAR")
}
if (!mask_kind %in% c("CCI", "GLC")) {
  stop_msg("Unsupported mask: ", mask_kind, ". Use CCI or GLC")
}

g_min <- cfg$luh2$pasture_mask$grass_min
p_min <- cfg$luh2$pasture_mask$pasture_min
beta <- cfg$luh2$pasture_mask$pasture_grass_ratio_min
luh_year_1 <- as.integer(cfg$luh2$pasture_mask$start_year)
luh_year_2 <- as.integer(cfg$luh2$pasture_mask$end_year)
nonveg_year <- as.integer(cfg$esa_cci$nonvegetated$year)
alpha_water <- cfg$esa_cci$nonvegetated$water_threshold
alpha_ice <- cfg$esa_cci$nonvegetated$ice_threshold
cci_persistence_years <- as.integer(cfg$esa_cci$used_land$persistence_years)
glc_persistence_years <- as.integer(cfg$glc$used_land$persistence_years)

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
  file.path(
    cfg$paths$masks_cci_dir,
    sprintf(
      "mask_used_frac_fused_alpha%s_k%d_%d-%d_0p05.tif",
      tok(cfg$esa_cci$used_land$threshold),
      cci_persistence_years,
      year_1,
      year_2
    )
  )
} else {
  file.path(
    cfg$paths$masks_glc_dir,
    sprintf(
      "mask_used_ge%d_%d-%d_0p05.tif",
      glc_persistence_years,
      year_1,
      year_2
    )
  )
}
if (!file.exists(mask_path)) {
  stop_msg("Missing land-use mask: ", mask_path)
}

drop_mask <- rast(mask_path)
drop_mask <- align_to_template(drop_mask, ref005, method = "near")
stopifnot(compareGeom(drop_mask, ref005, stopOnError = FALSE))

combine_or <- function(a, b) {
  # Combine drop masks; an NA component is treated as no additional exclusion.
  b <- align_to_template(b, a, method = "near")
  app(
    c(a, b),
    fun = function(v) {
      as.integer(any(v >= 1, na.rm = TRUE))
    }
  )
}

nonveg_dir <- file.path(cfg$paths$masks_root_dir, "mask_nonvegetated")
nonveg_path <- file.path(
  nonveg_dir,
  sprintf(
    "mask_nonvegetated_CCI_%d_alphaW%s_alphaI%s_0p05.tif",
    nonveg_year,
    tok(alpha_water),
    tok(alpha_ice)
  )
)
if (!file.exists(nonveg_path)) {
  stop_msg("Missing non-vegetated mask: ", nonveg_path)
}
nonveg <- rast(nonveg_path)
drop_mask <- combine_or(drop_mask, nonveg)


luh_dir <- file.path(cfg$paths$masks_root_dir, "mask_luh_overlap")
luh_name <- sprintf(
  "mask_luh_overlap_%s_Gmin%s_Pmin%s_beta%s_%d-%d_0p05_rep.tif",
  mask_kind,
  tok(g_min),
  tok(p_min),
  tok(beta),
  luh_year_1,
  luh_year_2
)
luh_path <- file.path(luh_dir, luh_name)
if (!file.exists(luh_path)) {
  stop_msg("Missing LUH2 pasture-overlap mask: ", luh_path)
}
luh <- rast(luh_path)
drop_mask <- combine_or(drop_mask, luh)
vals_ok <- try(all(values(drop_mask) %in% c(0, 1, NA)), silent = TRUE)
if (inherits(vals_ok, "try-error") || !isTRUE(vals_ok)) {
  stop_msg("Combined mask has values outside {0,1,NA}")
}

mask_combined_path <- file.path(out_dir, "combined_mask_0p05.tif")
mask_sources <- c(mask_path, nonveg_path, luh_path)
mask_needs_update <- !file.exists(mask_combined_path) ||
  any(file.mtime(mask_sources) > file.mtime(mask_combined_path))

if (mask_needs_update) {
  writeRaster(
    drop_mask,
    mask_combined_path,
    overwrite = TRUE,
    wopt = wopt_byte(FALSE, na = 255L)
  )
}
months <- format(
  seq(
    as.Date(sprintf("%d-01-01", cfg$project$years$lai_start)),
    as.Date(sprintf("%d-12-01", cfg$project$years$lai_end)),
    by = "month"
  ),
  "%Y%m"
)
files <- file.path(in_dir, sprintf("%s_%s_0p05.tif", var, months))
missing_files <- files[!file.exists(files)]
if (length(missing_files)) {
  stop("Missing expected monthly files:\n", paste(missing_files, collapse = "\n"))
}

wopt <- wopt_f32(FALSE)

for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_%s_0p05_masked.tif", var, ym))

  do_write <- !file.exists(out) ||
    file.mtime(f) > file.mtime(out) ||
    file.mtime(mask_combined_path) > file.mtime(out)

  if (do_write) {
    r <- rast(f)
    r_masked <- terra::mask(r,
      drop_mask,
      maskvalues = 1,
      updatevalue = NA
    )
    writeRaster(r_masked, out, overwrite = TRUE, wopt = wopt)
    rm(r, r_masked)
    gc(verbose = FALSE)
  }
}

gc()
