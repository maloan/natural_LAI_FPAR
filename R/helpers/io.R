# =============================================================================
# io.R — Helper functions for input/output operations, including GDAL options,
# raster reading/writing, and file path management.
# =============================================================================

gdal_co <- function(predictor = 2L,
                    speed_over_size = FALSE) {
  # Return GDAL creation options for GeoTIFF output.
  base <- c("TILED=YES", "BIGTIFF=IF_SAFER")
  comp <- if (isTRUE(speed_over_size)) {
    "COMPRESS=LZW"
  } else {
    "COMPRESS=DEFLATE"
  }
  c(
    base,
    comp,
    sprintf("PREDICTOR=%d", predictor),
    "NUM_THREADS=ALL_CPUS"
  )
}

gdal_co_int <- function(speed_over_size = FALSE) {
  # Return GDAL creation options for integer GeoTIFF output.
  gdal_co(predictor = 2L, speed_over_size = speed_over_size)
}

gdal_co_f32 <- function(speed_over_size = FALSE) {
  # Return GDAL creation options for float GeoTIFF output.
  gdal_co(predictor = 3L, speed_over_size = speed_over_size)
}

wopt_byte <- function(speed_over_size = FALSE,
                      na = 255L) {
  # Return a list of options for writing byte rasters with optional NA value.
  list(
    datatype = "INT1U",
    NAflag   = as.integer(na),
    gdal     = gdal_co_int(speed_over_size)
  )
}

wopt_int <- function(speed_over_size = FALSE,
                     na = NA_integer_) {
  # Return a list of options for writing integer rasters with optional NA value.
  out <- list(datatype = "INT2U", gdal = gdal_co_int(speed_over_size))
  if (!is.na(na)) {
    out$NAflag <- as.integer(na)
  }
  out
}

wopt_f32 <- function(speed_over_size = FALSE,
                     na = -9999) {
  #   Return a list of options for writing float rasters with optional NA value.
  list(
    datatype = "FLT4S",
    NAflag   = as.numeric(na),
    gdal     = gdal_co_f32(speed_over_size)
  )
}

read_trend <- function(path,
                       label = NULL,
                       template = NULL) {
  # Read a single-layer raster from the specified path, checking for existence
  # and optional template alignment.
  if (!file.exists(path)) {
    stop("Trend file not found: ", path, if (!is.null(label)) {
      paste0(" (", label, ")")
    })
  }
  r <- terra::rast(path)
  if (terra::nlyr(r) != 1) {
    stop(
      "Expected single-layer raster in ",
      path,
      if (!is.null(label)) {
        paste0(" (", label, ")")
      },
      ", got ",
      terra::nlyr(r)
    )
  }
  if (!is.null(template)) {
    terra::compareGeom(r, template, stopOnError = TRUE)
  }
  terra::values(r, dataframe = FALSE)
}

trend_path_factory <- function(var,
                               met,
                               source,
                               run_tag = NULL,
                               is_relative = FALSE) {
  # Return the file path for a trend raster based on the specified parameters.
  if (source == "unmasked") {
    suffix <- if (is_relative) {
      "trend_relative_peryear"
    } else {
      "trend_slope_peryear"
    }
    file.path(
      here::here("analysis", "unmasked", "0p25"),
      sprintf("%s_georef_%s_%s_0p25.nc", var, met, suffix)
    )
  } else {
    if (is.null(run_tag)) {
      stop("run_tag required for ", source, " source")
    }
    suffix <- if (is_relative) {
      "trend_relative_peryear"
    } else {
      "trend_slope_peryear"
    }
    file.path(
      here::here(
        "output",
        run_tag,
        "eval",
        sprintf("trend_%s_%s", var, source)
      ),
      sprintf("%s_%s_%s_0p25.nc", var, met, suffix)
    )
  }
}

analysis_raster_path <- function(var,
                                 met,
                                 source,
                                 run_tag = NULL,
                                 kind = c("metric", "trend_relative", "trend_slope", "trend_mk_pval"),
                                 grid_tag = "0p25") {
  # Return the file path for an analysis raster based on the specified parameters.
  kind <- match.arg(kind)

  if (source == "unmasked") {
    file_name <- switch(kind,
      metric = sprintf("%s_georef_%s_%s.nc", var, met, grid_tag),
      trend_relative = sprintf(
        "%s_georef_%s_trend_relative_peryear_%s.nc",
        var,
        met,
        grid_tag
      ),
      trend_slope = sprintf(
        "%s_georef_%s_trend_slope_peryear_%s.nc",
        var,
        met,
        grid_tag
      ),
      trend_mk_pval = sprintf("%s_georef_%s_trend_mk_pval_%s.nc", var, met, grid_tag)
    )

    return(file.path(here::here("analysis", "unmasked", grid_tag), file_name))
  }

  if (is.null(run_tag)) {
    stop("run_tag required for source ", source)
  }

  file_name <- switch(kind,
    metric = sprintf("%s_%s_%s.nc", var, met, grid_tag),
    trend_relative = sprintf("%s_%s_trend_relative_peryear_%s.nc", var, met, grid_tag),
    trend_slope = sprintf("%s_%s_trend_slope_peryear_%s.nc", var, met, grid_tag),
    trend_mk_pval = sprintf("%s_%s_trend_mk_pval_%s.nc", var, met, grid_tag)
  )

  file.path(here::here(
    "output",
    run_tag,
    "eval",
    sprintf("trend_%s_%s", var, source)
  ), file_name)
}


round_numeric <- function(df, digits = 5) {
  # Round all numeric columns in a data frame to the specified number of digits.
  dplyr::mutate(df, dplyr::across(dplyr::where(is.numeric), ~ round(.x, digits)))
}

load_checked_raster <- function(path,
                                template,
                                label = NULL,
                                n_layers = NULL,
                                first_layer = FALSE) {
  # Load a raster from the specified path, checking for existence, template
  # alignment, and optional layer count.
  if (!file.exists(path)) {
    stop("Missing raster", if (!is.null(label)) {
      paste0(" (", label, ")")
    }, ": ", path)
  }

  r <- terra::rast(path)
  terra::compareGeom(r[[1]], template, stopOnError = TRUE)

  if (!is.null(n_layers) && terra::nlyr(r) != n_layers) {
    stop(
      "Unexpected number of layers in ",
      path,
      ": expected ",
      n_layers,
      ", got ",
      terra::nlyr(r)
    )
  }

  if (isTRUE(first_layer)) {
    return(r[[1]])
  }

  r
}

safe_division <- function(num, den, positive_denominator = FALSE) {
  # Perform safe division, returning NA for non-finite results or if the
  # denominator is non-positive when specified.
  if (isTRUE(positive_denominator)) {
    return(ifelse(is.finite(den) & den > 0, num / den, NA_real_))
  }

  out <- num / den
  ifelse(is.finite(out), out, NA_real_)
}
tok <- function(x) {
  # Convert numeric values to a safe string representation for filenames.
  gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
}

extract_ym_from_filename <- function(p) {
  # Extract year-month string (YYYYMM) from a filename.
  s <- tools::file_path_sans_ext(basename(p))
  m <- regexpr("(19|20)\\d{2}(0[1-9]|1[0-2])", s, perl = TRUE)
  if (m[1] > 0) {
    substr(s, m[1], m[1] + attr(m, "match.length") - 1)
  } else {
    s
  }
}


filter_by_year_range <- function(values, year_min, year_max) {
  # Filter a vector of year values to those within the specified range.
  which(values >= year_min & values <= year_max)
}


align_to_template <- function(r, template, method = "bilinear") {
  # Align a raster to a template raster, resampling if necessary.
  if (!terra::compareGeom(r, template, stopOnError = FALSE)) {
    r <- terra::resample(r, template, method = method)
  }
  r
}

exp_ <- function(p) {
  # Expand and normalize a file path, handling user home directory and
  # converting to a standard format.
  normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
}

cfg_read <- function() {
  # Read configuration options from a YAML file specified by environment
  # variables, with fallbacks.
  root <- exp_(Sys.getenv("SNU_LAI_FPAR_ROOT", unset = "~/GitHub/natural_LAI_FPAR"))
  run_tag <- Sys.getenv("run_tag", "alpha_0.1")
  cfg_file <- file.path(root, "config", sprintf("config_%s.yml", run_tag))

  # Fall back to generic config if run_tag specific one doesn't exist
  if (!file.exists(cfg_file)) {
    cfg_file <- file.path(root, "config", "config.yml")
  }

  yaml::read_yaml(cfg_file)
}
