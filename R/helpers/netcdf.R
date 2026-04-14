## =============================================================================
## netcdf.R — NetCDF selection and raster alignment helpers
## =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

pick_varname <- function(nc, cfg_vars, var) {
  vars <- names(nc$var)
  primary <- cfg_vars$nc_var_name_primary %||% var
  fallback <- cfg_vars$nc_var_name_fallback %||% "auto_first_variable"

  var_u <- toupper(var)
  extras <- switch(var_u,
    LAI  = c("lai", "LAI", "Leaf_Area_Index"),
    FPAR = c("FPAR", "fPAR", "fpar", "Fpar"),
    character(0)
  )

  candidates <- unique(c(primary, fallback, extras))
  for (nm in candidates) {
    if (identical(nm, "auto_first_variable")) next
    if (nm %in% vars) {
      return(nm)
    }
  }

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

  vars[[1L]]
}

get_lonlat <- function(nc, vcfg) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for get_lonlat().", call. = FALSE)
  }

  loncand <- unique(c(vcfg$nc_lon_name, "lon", "longitude", "LONGITUDE"))
  latcand <- unique(c(vcfg$nc_lat_name, "lat", "latitude", "LATITUDE"))

  lonv <- latv <- NULL

  for (d in nc$dim) {
    dn <- tolower(d$name)
    if (is.null(lonv) && dn %in% tolower(loncand)) {
      lonv <- ncdf4::ncvar_get(nc, d$name)
    }
    if (is.null(latv) && dn %in% tolower(latcand)) {
      latv <- ncdf4::ncvar_get(nc, d$name)
    }
  }

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

align_to <- function(r, ref, method = c("near", "bilinear")) {
  method <- match.arg(method)
  if (terra::compareGeom(r, ref, stopOnError = FALSE)) {
    r
  } else {
    terra::resample(r, ref, method)
  }
}
