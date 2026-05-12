## =============================================================================
## bootstrap_ci.R — Standardized confidence interval calculation
##
## This module provides unified bootstrap CI functions used across analysis
## scripts. All functions use weighted resampling with replacement (Bayesian
## bootstrap) to account for spatial variation and area weighting in trend
## statistics.
##
## Methods:
##   - Global weighted mean CI: bootstrap_ci_global()
##   - By-class CI: bootstrap_ci_by_class()
##
## Statistical notes:
##   - Bootstrap uses area-weighted resampling with replacement (n = n)
##   - Confidence intervals are percentile-based (quantile type = 7)
##   - Fixed seed (42) ensures reproducibility across runs
##   - Default: n_boot = 1000, conf = 0.95 (adjust as needed)
##   - Significance: CI spans zero ⟹ not significant (* flag if CI > 0)
## =============================================================================

# Global weighted mean bootstrap confidence interval
# Args:
#   x: numeric vector of values (e.g., trend coefficients)
#   w: numeric vector of weights (e.g., area in km²)
#   n_boot: number of bootstrap samples
#   conf: confidence level (0.95 = 95% CI)
# Returns: list with mean, lower, upper, width, sig flag
bootstrap_ci_global <- function(x, w, n_boot = 1000L, conf = 0.95) {
  if (length(x) != length(w)) {
    stop("Length of x and w must match")
  }

  ok <- is.finite(x) & is.finite(w) & w > 0
  x_ok <- x[ok]
  w_ok <- w[ok]
  n <- length(x_ok)

  if (n < 2) {
    return(list(
      mean = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      width = NA_real_,
      sig = "",
      n_eff = n
    ))
  }

  # Normalize weights for probability resampling
  w_norm <- w_ok / sum(w_ok, na.rm = TRUE)

  # Bootstrap resampling with replacement
  boot_means <- numeric(n_boot)
  set.seed(42)

  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE, prob = w_norm)
    boot_means[i] <- sum(x_ok[idx] * w_ok[idx]) / sum(w_ok[idx])
  }

  alpha <- (1 - conf) / 2
  lower <- quantile(boot_means, alpha, na.rm = TRUE, type = 7)
  upper <- quantile(boot_means, 1 - alpha, na.rm = TRUE, type = 7)

  # Significance: CI does not cross zero
  is_sig <- lower * upper > 0

  list(
    mean = mean(boot_means, na.rm = TRUE),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    width = as.numeric(upper - lower),
    sig = ifelse(is_sig, "*", ""),
    n_eff = n
  )
}

# By-class bootstrap confidence intervals
# For stratified estimates (e.g., trends by land-cover class or climate zone)
# Args:
#   r_values: numeric vector of values per pixel (e.g., trend, zealmean, yearmax)
#   z_values: integer vector of zone/class IDs per pixel
#   w_values: numeric vector of weights per pixel (e.g., area)
#   n_boot: number of bootstrap samples
#   conf: confidence level
# Returns: tibble with one row per class
bootstrap_ci_by_class <- function(r_values, z_values, w_values,
                                  n_boot = 1000L, conf = 0.95) {
  if (length(r_values) != length(z_values) ||
    length(r_values) != length(w_values)) {
    stop("All input vectors must have same length")
  }

  classes <- sort(unique(z_values[!is.na(z_values)]))
  results <- list()

  set.seed(42)

  for (cls in classes) {
    idx <- which(
      z_values == cls &
        is.finite(r_values) &
        is.finite(w_values) &
        w_values > 0
    )

    if (length(idx) == 0) {
      results[[length(results) + 1]] <- tibble::tibble(
        class = cls,
        n_pixels = 0L,
        mean_est = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        ci_width = NA_real_,
        sig_flag = ""
      )
      next
    }

    r_cls <- r_values[idx]
    w_cls <- w_values[idx]
    n <- length(r_cls)

    # Normalized weights for probability resampling
    w_norm <- w_cls / sum(w_cls, na.rm = TRUE)

    # Point estimate
    wmean_est <- sum(r_cls * w_cls, na.rm = TRUE) / sum(w_cls, na.rm = TRUE)

    # Bootstrap resampling
    boot_means <- numeric(n_boot)
    for (i in seq_len(n_boot)) {
      boot_idx <- sample.int(n, size = n, replace = TRUE, prob = w_norm)
      boot_means[i] <- sum(r_cls[boot_idx] * w_cls[boot_idx], na.rm = TRUE) / sum(w_cls[boot_idx], na.rm = TRUE)
    }

    alpha <- (1 - conf) / 2
    ci_lower <- quantile(boot_means, alpha, na.rm = TRUE, type = 7)
    ci_upper <- quantile(boot_means, 1 - alpha, na.rm = TRUE, type = 7)

    # Significance check
    is_sig <- ci_lower * ci_upper > 0 & !is.na(ci_lower) & !is.na(ci_upper)

    results[[length(results) + 1]] <- tibble::tibble(
      class = cls,
      n_pixels = as.integer(n),
      mean_est = wmean_est,
      ci_lower = as.numeric(ci_lower),
      ci_upper = as.numeric(ci_upper),
      ci_width = as.numeric(ci_upper - ci_lower),
      sig_flag = ifelse(is_sig, "*", "")
    )
  }

  dplyr::bind_rows(results)
}
