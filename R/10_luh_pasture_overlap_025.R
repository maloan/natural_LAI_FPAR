## =============================================================================
# 11_apply_mask_0p05.R — Apply drop masks to monthly LAI/FPAR at 0.05°
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
# source(here("R","helpers", "geom.R"))
source(here("R", "helpers", "viz.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)

VAR <- toupper(Sys.getenv("VAR", "FPAR"))
MASK <- toupper(Sys.getenv("MASK", "CCI"))

SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
OVERWRITE <- as.logical(Sys.getenv("OVERWRITE", "FALSE"))
REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))
SRC <- toupper(Sys.getenv("GRASS_SOURCE", "CCI"))
GMIN <- as.numeric(Sys.getenv("G_MIN", "0.1"))
PMIN <- as.numeric(Sys.getenv("P_MIN", "0.1"))
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.50"))

ref005 <- rast(cfg$grids$grid_005$ref_raster)

in_dir <- if (VAR == "LAI") {
  cfg$paths$georef_lai_0p05_dir
} else {
  cfg$paths$georef_fpar_0p05_dir
}
out_dir <- switch(MASK,
  CCI = if (VAR == "LAI") {
    cfg$paths$masked_lai_cci_0p05_dir
  } else {
    cfg$paths$masked_fpar_cci_0p05_dir
  },
  GLC = if (VAR == "LAI") {
    cfg$paths$masked_lai_glc_0p05_dir
  } else {
    cfg$paths$masked_fpar_glc_0p05_dir
  }
)

stopifnot(is.character(out_dir), length(out_dir) == 1, nzchar(out_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ql_dir <- file.path(out_dir, "quicklooks")
dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

Y1 <- cfg$project$years$cci_start
Y2 <- cfg$project$years$cci_end

mask_path <- if (MASK == "CCI") {
  find_one(
    cfg$paths$masks_cci_dir,
    sprintf("mask_used_.*_%d-%d_0p05\\.tif$", Y1, Y2)
  )
} else {
  N <- as.integer(Sys.getenv("USED_N_YEARS", "3"))
  find_one(
    cfg$paths$masks_glc_dir,
    sprintf("mask_used_ge%d_%d-%d_0p05\\.tif$", N, Y1, Y2)
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
  abi <- rast(abi_path[order(file.info(abi_path)$mtime,
    decreasing = TRUE
  )][1])
  mask <- combine_or(mask, abi)
}

rx <- sprintf(
  "mask_luh_overlap_%s_Gmin%s_Pmin%s_alpha%s_%d-%d_0p05_rep\\.tif$",
  SRC, tok(GMIN), tok(PMIN), tok(ALPHA), Y1, Y2
)
luh_path <- find_one(
  file.path(cfg$paths$masks_root_dir, "mask_luh_overlap"), rx
)
luh <- rast(luh_path)
mask <- combine_or(mask, luh)
vals_ok <- try(all(values(mask) %in% c(0, 1, NA)), silent = TRUE)

mask_combined_path <- file.path(out_dir, "combined_mask_0p05.tif")
if (!file.exists(mask_combined_path) || OVERWRITE) {
  writeRaster(mask, mask_combined_path,
    overwrite = TRUE,
    wopt = wopt_byte(opts$SPEED_OVER_SIZE, na = 255L)$wopt
  )
}

files <- sort(list.files(
  in_dir,
  pattern = paste0("^", VAR, "_\\d{6}_0p05\\.tif$"), full.names = TRUE
))
stopifnot(length(files) > 0L)

wopt <- wopt_f32(opts$SPEED_OVER_SIZE)$wopt

for (f in files) {
  ym <- extract_ym_from_filename(f)
  out <- file.path(out_dir, sprintf("%s_%s_0p05_masked.tif", VAR, ym))
  ql1 <- file.path(ql_dir, sprintf("quicklook_%s_before_after_%s.png", VAR, ym))
  ql2 <- file.path(ql_dir, sprintf("quicklook_%s_after_%s.png", VAR, ym))

  do_write <- OVERWRITE || !file.exists(out) || !SKIP_EXISTING
  do_ql <- substr(ym, 5, 6) %in% c("01", "07") &&
    (REMAKE_QL || !file.exists(ql1) || !file.exists(ql2))

  if (!do_write && !do_ql) {
    next
  }

  r <- rast(f)

  if (do_write) {
    r_masked <- mask(r, mask, maskvalues = 1, updatevalue = NA)
    writeRaster(r_masked, out, overwrite = TRUE, wopt = wopt)
  } else {
    r_masked <- rast(out)
  }

  if (do_ql) {
    if (REMAKE_QL || !file.exists(ql1)) {
      quicklook_before_after(r, r_masked, ym, title = VAR, ql_dir = ql_dir)
    }
    if (REMAKE_QL || !file.exists(ql2)) {
      quicklook_after_full(r_masked, ym, title = VAR, ql_dir = ql_dir)
    }
  }
}

gc()
message("Masked ", VAR, " written to: ", out_dir)
