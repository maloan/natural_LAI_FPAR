## =============================================================================
# 01_georef_0p05.R — Georeference monthly LAI/FPAR to 0.05° grid
## =============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(here)
})

# --- repo helpers ---
ROOT <- here()

source(here("R", "utils.R"))
source(here("R", "io.R"))
source(here("R", "geom.R"))
source(here("R", "viz.R"))
source(here("R", "options.R"))

cfg  <- cfg_read()
opts <- opts_read()

terraOptions(progress = 1, memfrac = 0.25)

REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))
wopt <- wopt_f32(opts$SPEED_OVER_SIZE)

# reference grid
ref005 <- rast(cfg$grids$grid_005$ref_raster)

# --- choose variable ---
VAR <- toupper(Sys.getenv("VAR", "LAI"))
stopifnot(VAR %in% c("LAI", "FPAR"))
Vcfg <- cfg$variables[[tolower(VAR)]]

# plotting
pal_green <- function(n = 64) hcl.colors(n, palette = "Greens", rev = TRUE)
col_na <- "#bdbdbd"

# --- paths from config ---
if (VAR == "LAI") {
  out_georef <- cfg$paths$georef_lai_0p05_dir
  out_quick  <- file.path(out_georef, "quicklooks")
  nc_dir     <- cfg$paths$lai_nc_dir
  zlim_plot  <- c(0, 1)
} else {
  out_georef <- cfg$paths$georef_fpar_0p05_dir
  out_quick  <- file.path(out_georef, "quicklooks")
  nc_dir     <- cfg$paths$fpar_nc_dir
  zlim_plot  <- c(0, 1)
}

dir.create(out_georef, TRUE, showWarnings = FALSE)
dir.create(out_quick, TRUE, showWarnings = FALSE)

# --- discover inputs ---
stopifnot(dir.exists(nc_dir))
files <- sort(list.files(nc_dir, pattern = "\\.nc(4)?$", full.names = TRUE))
stopifnot(length(files) > 0L)

# --- helpers unchanged ---------------------------------------------------------
pick_varname <- function(nc, cfg_vars, VAR) {
  vars <- names(nc$var)
  primary  <- cfg_vars$nc_var_name_primary %||% VAR
  fallback <- cfg_vars$nc_var_name_fallback %||% "auto_first_variable"
  candidates <- if (VAR == "LAI") {
    unique(c(primary, fallback, "lai", "LAI", "Leaf_Area_Index"))
  } else {
    unique(c(primary, fallback, "FPAR", "fPAR", "fpar", "Fpar"))
  }
  for (nm in candidates) if (nm %in% vars) return(nm)
  vars[[1L]]
}

extract_ym_from_filename <- function(p) {
  s <- tools::file_path_sans_ext(basename(p))
  m <- regexpr("(19|20)\\d{2}(0[1-9]|1[0-2])", s, perl = TRUE)
  if (m[1] > 0) substr(s, m[1], m[1] + attr(m, "match.length") - 1) else s
}

# --- processing loop -----------------------------------------------------------
for (f in files) {

  ym <- extract_ym_from_filename(f)

  out_tif <- file.path(out_georef, sprintf("%s_%s_0p05.tif", VAR, ym))
  out_png <- file.path(out_quick,  sprintf("%s_%s_0p05.png", VAR, ym))

  if (REMAKE_QL && file.exists(out_tif)) {
    r <- rast(out_tif)
    png(out_png, 1100, 550, res = 96)
    par(mar = c(3,3,2,6))
    plot(
      r, main = sprintf("%s %s (0.05°)", VAR, ym),
      zlim = zlim_plot, col = pal_green(64), colNA = col_na, axes = TRUE
    )
    dev.off()
    next
  }

  if (!opts$FORCE && opts$SKIP_EXISTING && file.exists(out_tif))
    next

  nc <- ncdf4::nc_open(f)
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE), add = TRUE)

  varnm <- pick_varname(nc, Vcfg, VAR)

  arr <- ncdf4::ncvar_get(nc, varnm)
  arr <- extract_time_slice(arr, varnm, nc, ym)
  if (is.null(arr))
    stop(sprintf("Time slice '%s' missing in '%s' (%s).", ym, basename(f), varnm))

  # lon/lat discovery
  loncand <- c("lon","longitude","LONGITUDE")
  latcand <- c("lat","latitude","LATITUDE")
  lonv <- latv <- NULL
  for (d in nc$dim) {
    dn <- tolower(d$name)
    if (dn %in% tolower(loncand)) lonv <- ncdf4::ncvar_get(nc, d$name)
    if (dn %in% tolower(latcand)) latv <- ncdf4::ncvar_get(nc, d$name)
  }
  if (length(dim(arr)) == 2L && !is.null(lonv) && !is.null(latv))
    arr <- maybe_transpose_lonlat(arr, length(lonv), length(latv))

  # attributes
  atts <- ncdf4::ncatt_get(nc, varnm)
  fillv <- if (!is.null(atts$`_FillValue`)) atts$`_FillValue`
  else if (!is.null(atts$missing_value)) atts$missing_value else NA
  if (!is.na(fillv)) arr[arr == fillv] <- NA_real_

  if (!is.null(atts$scale_factor) && !identical(atts$scale_factor, 1)) {
    rng <- range(arr, na.rm = TRUE)
    if (rng[2] > 1.5) arr <- arr * as.numeric(atts$scale_factor)
  }
  if (!is.null(atts$add_offset) && !identical(atts$add_offset, 0))
    arr <- arr + as.numeric(atts$add_offset)

  # clamp + auto FPAR scaling
  if (!is.null(Vcfg$clamp)) {
    lo <- Vcfg$clamp$min
    hi <- Vcfg$clamp$max
    if (VAR == "FPAR") {
      rng <- range(arr[is.finite(arr)], na.rm = TRUE)
      if (rng[2] > 1.5 && hi <= 1) arr <- arr / 100
    }
    arr[arr < lo] <- lo
    arr[arr > hi] <- hi
  }

  # raster → georeg
  r <- rast(arr)
  ext(r) <- ext(-180,180,-90,90)
  crs(r) <- crs(ref005)
  r <- rotate_if_360(r, lon_vals = lonv)

  if (!same_grid(r, ref005)) {
    message(sprintf("Resampling %s → 0.05° grid (bilinear)", basename(f)))
    r <- resample(r, ref005, method = "bilinear")
  }

  writeRaster(r, out_tif, overwrite = TRUE, wopt = wopt)

  # plot
  if (REMAKE_QL || !file.exists(out_png)) {
    png(out_png, 1100, 550, res = 96)
    par(mar = c(3,3,2,6))
    plot(
      r,
      main = sprintf("%s %s (0.05°)", VAR, ym),
      zlim = zlim_plot,
      col = pal_green(64),
      colNA = col_na,
      axes = TRUE
    )
    dev.off()
  }
}

gc()
message("Done 01_georef_monthlies_to_0p05.R")
