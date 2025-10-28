## =============================================================================
# 01_georef_0p05.R — Georeference monthly LAI/FPAR to 0.05° grid
#
# Purpose
#   Read monthly LAI or FPAR from NetCDF files, apply CF scaling and clamping,
#   orient and rotate to lon–lat coordinates, resample to the global 0.05°
#   reference grid, and export as GeoTIFFs with optional quicklook PNGs.
#
# Inputs
#   - NetCDF files: cfg$paths$lai_nc_dir or cfg$paths$fpar_nc_dir
#   - Reference grid: cfg$grids$grid_005$ref_raster (EPSG:4326)
#   - Config: cfg_read(), opts_read()
#
# Outputs
#   - GeoTIFFs: {VAR}_{YYYYMM}_0p05.tif in georef_{VAR}_0p05_dir
#   - Quicklooks: {VAR}_{YYYYMM}_0p05.png in quicklooks/
#
# Environment variables
#   VAR        (string, default "FPAR") — variable to process ("LAI" or "FPAR")
#   REMAKE_QL  (logical, default FALSE) — force remake of quicklooks
#
# Dependencies
#   Packages: ncdf4, terra
#   Sourced helpers: utils.R, io.R, geom.R, viz.R, options.R
#
# Processing overview
#   1) Discover monthly NetCDF files.
#   2) Read variable and time slice, apply scale/offset/fill corrections.
#   3) Clamp and rescale values; handle 0–360° longitude.
#   4) Snap to 0.05° reference grid (bilinear).
#   5) Write GeoTIFF and optional PNG quicklook.
## =============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
})

# --- repo helpers ---
# set root
ROOT <- normalizePath(path.expand(
  Sys.getenv("SNU_LAI_FPAR_ROOT", "~/GitHub/natural_LAI_FPAR")
),
winslash = "/",
mustWork = FALSE)
# source other files
source(file.path(ROOT, "R", "utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))
source(file.path(ROOT, "R", "viz.R"))
source(file.path(ROOT, "R", "options.R"))
cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)
# list of GDAL/terra options for a 32-bit float GeoTIFF (datatype, tiling,
# compression, etc.)
REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))
wopt <- wopt_f32(opts$SPEED_OVER_SIZE)

# define raster
ref005 <- rast(path.expand(cfg$grids$grid_005$ref_raster))
# --- choose variable: VAR=LAI|FPAR (env) ---
VAR <- toupper(Sys.getenv("VAR", "FPAR"))
stopifnot(VAR %in% c("LAI", "FPAR"))
Vcfg <- cfg$variables[[tolower(VAR)]]

# ---- plotting palette (sequential green: low=light, high=dark; NA=grey) ----
pal_green <- function(n = 64)
  hcl.colors(n, palette = "Greens", rev = TRUE)
col_na <- "#bdbdbd"  # grey for NA

# --- paths from config keyed by VAR ---
if (VAR == "LAI") {
  # Set directories and plotting setting
  out_georef <- path.expand(cfg$paths$georef_lai_0p05_dir)
  out_quick  <- file.path(path.expand(out_georef), "quicklooks")
  nc_dir     <- path.expand(cfg$paths$lai_nc_dir)
  zlim_plot  <- c(0, 8)
} else {
  # Set directories and plotting setting
  out_georef <- path.expand(cfg$paths$georef_fpar_0p05_dir)
  out_quick  <- file.path(path.expand(out_georef), "quicklooks")
  nc_dir     <- path.expand(cfg$paths$fpar_nc_dir)
  zlim_plot  <- c(0, 1)
}
# create respective directories
dir.create(out_georef, TRUE, showWarnings = FALSE)
dir.create(out_quick, TRUE, showWarnings = FALSE)

# --- discover inputs ---
stopifnot(dir.exists(nc_dir))
files <- sort(list.files(nc_dir, pattern = "\\.nc(4)?$", full.names = TRUE))
stopifnot(length(files) > 0L)

# --- options ---
pick_varname <- function(nc, cfg_vars, VAR) {
  #' Pick NetCDF variable name
  #'
  #' Selects the variable to read from an open NetCDF file by checking
  #' config-specified names and common aliases for `"LAI"` or `"FPAR"`. Returns
  #' the first match found, or the first variable in the file if none match.
  #'
  #' @param nc NetCDF handle from [ncdf4::nc_open()].
  #' @param cfg_vars Config list, may include `nc_var_name_primary` and
  #'   `nc_var_name_fallback`.
  #' @param VAR Character: `"LAI"` or `"FPAR"`.
  #' @return Character string with the chosen variable name.
  vars <- names(nc$var)
  primary  <- cfg_vars$nc_var_name_primary %||% VAR
  fallback <- cfg_vars$nc_var_name_fallback %||% "auto_first_variable"
  candidates <- if (VAR == "LAI") {
    unique(c(primary, fallback, "lai", "LAI", "Leaf_Area_Index"))
  } else {
    unique(c(primary, fallback, "FPAR", "fPAR", "fpar", "Fpar"))
  }
  for (nm in candidates)
    if (nm %in% vars) {
      return(nm)
    }
  vars[[1L]]
}

extract_ym_from_filename <- function(p) {
  s <- tools::file_path_sans_ext(basename(p))
  m <- regexpr("(19|20)\\d{2}(0[1-9]|1[0-2])", s, perl = TRUE)
  if (m[1] > 0) {
    substr(s, m[1], m[1] + attr(m, "match.length") - 1)
  } else {
    s
  }
}

# --- processing ---
for (f in files) {
  # extract year from filename
  ym <- extract_ym_from_filename(f)
  # construct output name
  out_tif <- file.path(out_georef, sprintf("%s_%s_0p05.tif", VAR, ym))
  out_png <- file.path(out_quick, sprintf("%s_%s_0p05.png", VAR, ym))

  if (REMAKE_QL && file.exists(out_tif)) {
    # Remake quicklooks and output tif exists
    # only remake plots
    r <- rast(out_tif)  # reuse georeffed output
    png(out_png, 1100, 550, res = 96)
    par(mar = c(3, 3, 2, 6))
    plot(
      r,
      main  = sprintf("%s %s (0.05°)", VAR, ym),
      zlim  = zlim_plot,
      col   = pal_green(64),
      colNA = col_na,
      axes  = TRUE
    )
    dev.off()
    next
  }
  if (!opts$FORCE &&
      opts$SKIP_EXISTING && file.exists(out_tif)) {
    # Do nothing
    next
  }

  nc <- ncdf4::nc_open(f) # Open nc file
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE)
          , add = TRUE) # close on exit

  # Pick variable name
  varnm <- pick_varname(nc, Vcfg, VAR)

  # array + time slice
  # reads full NetCDF for varname into R
  arr <- ncdf4::ncvar_get(nc, varnm)
  # subsets array to a single year, using the file’s time coordinate/metadata.
  # Returns the 2D (lon × lat) slice
  arr <- extract_time_slice(arr, varnm, nc, ym)
  if (is.null(arr)) {
    # If arraz is empty, stops
    stop(sprintf(
      "Time slice '%s' not found in '%s' for '%s'.",
      ym,
      basename(f),
      varnm
    ))
  }

  #lon/lat (optional) + orientation
  # find right variable name
  loncand <- c("lon", "longitude", "LONGITUDE")
  latcand <- c("lat", "latitude", "LATITUDE")
  # define lonv and latv as Null
  lonv <- latv <- NULL
  # find longitude and latitude dimensions
  for (d in nc$dim) {
    # For each match, read coordinate vector and store them in lonv and latv
    dn <- tolower(d$name)
    if (dn %in% tolower(loncand))
    {
      lonv <- ncdf4::ncvar_get(nc, d$name)
    }
    if (dn %in% tolower(latcand))
    {
      latv <- ncdf4::ncvar_get(nc, d$name)
    }
  }
  if (length(dim(arr)) == 2L && !is.null(lonv) && !is.null(latv)) {
    # Check the dimension siyes against length of lonv and latv, transpose if
    # needed for first dim to match lon and the second to match lat.
    arr <- maybe_transpose_lonlat(arr, lon_len = length(lonv), lat_len = length(latv))
  }

  # scale/offset/fill
  # read all attributes for varnm
  atts  <- ncdf4::ncatt_get(nc, varnm)
  # determine the fill value
  fillv <- if (!is.null(atts$`_FillValue`)) {
    atts$`_FillValue`
  } else if (!is.null(atts$missing_value)) {
    atts$missing_value
  } else {
    NA
  }
  if (!is.na(fillv)) {
    # replace fill codes in arr with NA_real_
    arr[arr == fillv] <- NA_real_
  }
  # Apply CF-style scaling
  if (!is.null(atts$scale_factor) &&
      !identical(atts$scale_factor, 1)) {
    # if scale_factor ≠ 1, multiply
    arr <- arr * as.numeric(atts$scale_factor)
  }
  if (!is.null(atts$add_offset)   &&
      !identical(atts$add_offset, 0)) {
    # if add_offset ≠ 0, add
    arr <- arr + as.numeric(atts$add_offset)
  }

  # clamp / auto-scale FPAR (0–100 → 0–1)
  if (!is.null(Vcfg$clamp)) {
    # if clamp range defines in config, read that
    lo <- as.numeric(Vcfg$clamp$min)
    hi <- as.numeric(Vcfg$clamp$max)
    if (VAR == "FPAR") {
      # if the data look like 0–100 (max > 1.5) but the configured range is ≤ 1,
      # divide by 100 to convert to 0–1
      rng <- range(arr[is.finite(arr)], na.rm = TRUE)
      if (rng[2] > 1.5 && hi <= 1)
        arr <- arr / 100
    }
    # clip values to lo and hi
    arr[arr < lo] <- lo
    arr[arr > hi] <- hi
  }

  # rasterize + georeg
  # make a terra raster from the 2-D array
  r <- rast(arr)
  # set geographic extent to full world lon/lat
  ext(r) <- ext(-180, 180, -90, 90)
  # assign the same coordinate reference system as ref005 (WGS84 lon/lat)
  crs(r) <- crs(ref005)
  # if longitudes are 0–360, rotate to −180–180 (dateline-safe), using lonv to
  # detect
  r <- rotate_if_360(r, lon_vals = lonv)

  # snap to 0.05° grid
  if (!same_grid(r, ref005)) {
    # if extent, resolution, CRS, or alignment differ interpolate r onto
    # ref005’s 0.05° grid using bilinear interpolation
    message(sprintf("Resampling %s → 0.05° grid (bilinear)", basename(f)))
    r <- resample(r, ref005, method = "bilinear")
  }

  # write
  writeRaster(r, out_tif, overwrite = TRUE, wopt = wopt)

  # quicklook
  if (REMAKE_QL || !file.exists(out_png)) {
    png(out_png, 1100, 550, res = 96)
    par(mar = c(3, 3, 2, 6))
    plot(
      r,
      main  = sprintf("%s %s (0.05°)", VAR, ym),
      zlim  = zlim_plot,
      col   = pal_green(64),
      colNA = col_na,
      # <- grey for NA
      axes  = TRUE
    )
    dev.off()
  }

}
gc()
message("Done 01_georef_monthlies_to_0p05.R")
