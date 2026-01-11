## =============================================================================
## utils.R — General-purpose helpers (config, tokens, IO discovery)
## =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------------------
# Paths + config
# ------------------------------------------------------------------------------

exp_ <- function(p) {
  normalizePath(path.expand(p),
    winslash = "/",
    mustWork = FALSE
  )
}

cfg_read <- function() {
  ROOT <- exp_(Sys.getenv(
    "SNU_LAI_FPAR_ROOT",
    unset = "~/GitHub/natural_LAI_FPAR"
  ))
  yaml::read_yaml(file.path(ROOT, "config", "config.yml"))
}

# ------------------------------------------------------------------------------
# Small string / filename helpers
# ------------------------------------------------------------------------------

tok <- function(x) gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))

extract_ym_from_filename <- function(p) {
  s <- tools::file_path_sans_ext(basename(p))
  m <- regexpr("(19|20)\\d{2}(0[1-9]|1[0-2])", s, perl = TRUE)
  if (m[1] > 0) {
    substr(s, m[1], m[1] + attr(m, "match.length") - 1)
  } else {
    s
  }
}
pick_varname <- function(nc, cfg_vars, VAR) {
  vars <- names(nc$var)
  primary <- cfg_vars$nc_var_name_primary %||% VAR
  fallback <- cfg_vars$nc_var_name_fallback %||% "auto_first_variable"

  var_u <- toupper(VAR)
  extras <- switch(var_u,
    LAI  = c("lai", "LAI", "Leaf_Area_Index"),
    FPAR = c("FPAR", "fPAR", "fpar", "Fpar"),
    character(0)
  )

  candidates <- unique(c(primary, fallback, extras))
  # 1) try explicit candidates (except the sentinel)
  for (nm in candidates) {
    if (identical(nm, "auto_first_variable")) next
    if (nm %in% vars) {
      return(nm)
    }
  }

  # 2) robust fallback: first non-coordinate var with >=2 dims
  coord_names <- tolower(c(
    "lon", "longitude", "x", "lat", "latitude",
    "y", "time", "time_counter", "t", "month"
  ))
  ok <- Filter(function(nm) {
    v <- nc$var[[nm]]
    if (is.null(v)) {
      return(FALSE)
    }
    if (tolower(nm) %in% coord_names) {
      return(FALSE)
    }
    nd <- length(v$dim)
    nd >= 2
  }, vars)

  if (length(ok)) {
    return(ok[[1]])
  }

  # 3) last resort: keep old behaviour
  vars[[1L]]
}


find_one <- function(dir, pattern) {
  cand <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(cand)) {
    stop("No file matching pattern '", pattern, "' in: ", dir, call. = FALSE)
  }
  cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]
}

# # ------------------------------------------------------------------------------
# # NetCDF lon/lat extraction (requires ncdf4 handle)
# # ------------------------------------------------------------------------------
#
get_lonlat <- function(nc, Vcfg) {
  # nc is expected to be an ncdf4 object (already opened)
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for get_lonlat().", call. = FALSE)
  }

  loncand <- unique(c(Vcfg$nc_lon_name, "lon", "longitude", "LONGITUDE"))
  latcand <- unique(c(Vcfg$nc_lat_name, "lat", "latitude", "LATITUDE"))

  lonv <- latv <- NULL

  # Try dimensions first
  for (d in nc$dim) {
    dn <- tolower(d$name)
    if (is.null(lonv) && dn %in% tolower(loncand)) {
      lonv <- ncdf4::ncvar_get(nc, d$name)
    }
    if (is.null(latv) && dn %in% tolower(latcand)) {
      latv <- ncdf4::ncvar_get(nc, d$name)
    }
  }

  # Fall back to variables
  if (is.null(lonv)) {
    for (nm in loncand) {
      if (!is.null(nc$var[[nm]])) {
        lonv <- ncdf4::ncvar_get(nc, nm)
        break
      }
    }
  }
  if (is.null(latv)) {
    for (nm in latcand) {
      if (!is.null(nc$var[[nm]])) {
        latv <- ncdf4::ncvar_get(nc, nm)
        break
      }
    }
  }
  # Fall back to variables
  if (is.null(lonv)) {
    for (nm in loncand) {
      if (!is.null(nc$var[[nm]])) {
        lonv <- ncdf4::ncvar_get(nc, nm)
        break
      }
    }
  }
  if (is.null(latv)) {
    for (nm in latcand) {
      if (!is.null(nc$var[[nm]])) {
        latv <- ncdf4::ncvar_get(nc, nm)
        break
      }
    }
  }

  list(lon = lonv, lat = latv)
}

# ------------------------------------------------------------------------------
# Raster helpers
# ------------------------------------------------------------------------------
gdal_wopt <- function(dtype = "FLT4S",
                      compress = "DEFLATE",
                      threads = "AUTO") {
  thr <- if (identical(threads, "AUTO")) {
    "NUM_THREADS=ALL_CPUS"
  } else {
    sprintf("NUM_THREADS=%s", threads)
  }
}

align_to <- function(r, ref, method = c("near", "bilinear")) {
  method <- match.arg(method)
  if (terra::compareGeom(r, ref, stopOnError = FALSE)) {
    r
  } else {
    terra::resample(r, ref, method)
  }
}
