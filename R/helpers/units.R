## =============================================================================
## units.R — Standardized unit conversion and validation
##
## Centralized location for all unit conversions (% per year, absolute per year)
## and unit validation across analysis scripts.
## =============================================================================

# Unit conversion: fractional to percentage
frac_to_pct <- function(x, scale_factor = 100) {
  # Convert from fractional (0-1) to percentage
  # Args: x = numeric vector, scale_factor = conversion (default 100 for %)
  x * scale_factor
}

# Unit conversion: percentage to fractional
pct_to_frac <- function(x) {
  # Convert from percentage to fractional (0-1)
  x / 100
}

# Get unit string for a variable
get_unit_label <- function(var, trend_type = "absolute") {
  # Args:
  #   var: "LAI" or "FPAR"
  #   trend_type: "absolute" or "relative"
  # Returns: character string for plot labels

  stopifnot(var %in% c("LAI", "FPAR"))
  stopifnot(trend_type %in% c("absolute", "relative"))

  if (trend_type == "absolute") {
    # Absolute trends in native units per year
    if (var == "LAI") {
      "LAI units yr⁻¹"
    } else {
      "fAPAR units yr⁻¹"
    }
  } else {
    # Relative trends as percentage per year
    "% yr⁻¹"
  }
}

# Validate unit consistency in data
validate_units <- function(x, var, trend_type = "absolute",
                           name = "trend values") {
  # Check that trend values are in expected range for given variable and type
  # Args:
  #   x: numeric vector of trend values
  #   var: "LAI" or "FPAR"
  #   trend_type: "absolute" or "relative"
  #   name: description for error messages
  # Returns: invisible(x) if valid, stop() if invalid

  stopifnot(var %in% c("LAI", "FPAR"))
  stopifnot(trend_type %in% c("absolute", "relative"))

  ok <- is.finite(x)
  if (sum(ok) == 0) {
    return(invisible(x)) # All NA, no check needed
  }

  x_ok <- x[ok]

  # Define expected ranges
  expected_ranges <- list(
    LAI = list(absolute = c(-0.1, 0.1), relative = c(-0.5, 0.5)),
    FPAR = list(absolute = c(-0.05, 0.05), relative = c(-0.5, 0.5))
  )

  range <- expected_ranges[[var]][[trend_type]]

  # Check extremes (1st and 99th percentiles)
  p01 <- quantile(x_ok, 0.01, na.rm = TRUE, type = 7)
  p99 <- quantile(x_ok, 0.99, na.rm = TRUE, type = 7)

  # Allow 50% beyond expected (sometimes data has outliers)
  if (p01 < range[1] * 1.5 || p99 > range[2] * 1.5) {
    warning(
      sprintf(
        "Unit validation warning for %s (%s, %s trend):\n",
        name, var, trend_type
      ),
      sprintf("  1st-99th percentile: [%.4f, %.4f]\n", p01, p99),
      sprintf("  Expected range: [%.4f, %.4f]\n", range[1], range[2]),
      "  Consider verifying unit conversion was applied correctly."
    )
  }

  invisible(x)
}

# Format trend values for CSV output (with units in header)
format_trend_csv <- function(df, var, trend_type = "absolute",
                             n_decimal = 4) {
  # Formats trend columns with consistent decimal places
  # Args:
  #   df: data frame with trend columns
  #   var: "LAI" or "FPAR" (for unit validation)
  #   trend_type: "absolute" or "relative"
  #   n_decimal: number of decimal places to round
  # Returns: data frame with formatted columns

  # Identify trend columns (common patterns)
  trend_cols <- grep("trend|mean|slope", names(df), ignore.case = TRUE)

  if (length(trend_cols) == 0) {
    warning("No trend columns identified in data frame")
    return(df)
  }

  # Round trend columns
  df[trend_cols] <- lapply(df[trend_cols], function(x) {
    if (is.numeric(x)) round(x, n_decimal) else x
  })

  invisible(df)
}

# Create CSV metadata comment header
make_csv_metadata <- function(variable, trend_type = "absolute",
                              mask = NULL, tau = NULL,
                              temporal_range = c(1982, 2024),
                              method = "Bayesian bootstrap",
                              n_boot = 1000, conf = 0.95) {
  # Generate standardized metadata header for CSV files
  # Args:
  #   variable: "LAI" or "FPAR"
  #   trend_type: "absolute" or "relative"
  #   mask: mask source (e.g., "CCI", "GLC", or NULL for unmasked)
  #   tau: tau parameter (e.g., "tau_0.1")
  #   temporal_range: c(year_start, year_end)
  #   method: "Bayesian bootstrap" or other CI method
  #   n_boot: number of bootstrap samples
  #   conf: confidence level (0.95 = 95%)
  # Returns: character vector of comment lines

  lines <- character()

  lines[1] <- "# METADATA: Global/regional trend summary"
  lines[2] <- sprintf("# Variable: %s", variable)
  lines[3] <- sprintf("# Trend type: %s (%s)", trend_type, get_unit_label(variable, trend_type))

  if (!is.null(mask)) {
    lines[4] <- sprintf("# Mask: %s", mask)
  }
  if (!is.null(tau)) {
    lines[5] <- sprintf("# Run tag (tau): %s", tau)
  }

  lines[6] <- sprintf(
    "# Temporal coverage: %d-%d (%d years)",
    temporal_range[1], temporal_range[2],
    temporal_range[2] - temporal_range[1] + 1
  )

  lines[7] <- sprintf(
    "# CI method: %s (n_boot=%d, confidence=%d%%)",
    method, n_boot, round(conf * 100)
  )

  lines[8] <- "# Significance: * = CI does not cross zero, blank = not significant"
  lines[9] <- "# Sample size: n_pixels = number of valid pixels per aggregation"
  lines[10] <- "# Area: area_km2 = total area in km² within valid domain"

  paste(lines, collapse = "\n")
}
