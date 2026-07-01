# =============================================================================
# bootstrap_ci.R — Spatial block bootstrap confidence intervals for weighted
# means
# =============================================================================

library(dplyr)

make_block_id <- function(area_raster, block_size_deg = 5) {
  # Create a spatial block ID for each pixel in the raster based on its
  # geographic coordinates.
  xy <- terra::crds(area_raster, df = TRUE, na.rm = FALSE)
  block_lon <- floor((xy[, 1] + 180) / block_size_deg)
  block_lat <- floor((xy[, 2] + 90) / block_size_deg)
  as.character(paste(block_lon, block_lat, sep = "_"))
}
# -------------------------------------------------------------------------
# Spatial block bootstrap for one weighted mean
# -------------------------------------------------------------------------
bootstrap_ci_global <- function(x,
                                w,
                                block_id,
                                n_boot = 1000L,
                                conf = 0.95) {
  # Compute a spatial block bootstrap confidence interval for a weighted mean of
  # x with weights w, using block_id to define spatial blocks.
  if (!(length(x) == length(w) &&
    length(x) == length(block_id))) {
    stop("x, w and block_id must have equal length.")
  }

  # Keep only valid weighted observations in valid spatial blocks.
  ok <- is.finite(x) & is.finite(w) & w > 0 & !is.na(block_id)
  dat <- data.frame(x = x[ok], w = w[ok], block = block_id[ok])
  if (nrow(dat) == 0) {
    return(
      list(
        mean  = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        width = NA_real_,
        sig   = "",
        n_eff = 0L
      )
    )
  }
  blocks <- unique(dat$block)
  n_blocks <- length(blocks)
  if (n_blocks < 2) {
    return(
      list(
        mean  = weighted.mean(dat$x, dat$w),
        lower = NA_real_,
        upper = NA_real_,
        width = NA_real_,
        sig   = "",
        n_eff = n_blocks
      )
    )
  }
  # Build an index of rows per block for efficient resampling.
  make_block_index <- function(block_id) {
    split(seq_along(block_id), block_id)
  }
  set.seed(26)
  block_idx <- make_block_index(dat$block)
  boot_est <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    sampled_blocks <- sample(names(block_idx), n_blocks, replace = TRUE)
    idx <- unlist(block_idx[sampled_blocks], use.names = FALSE)
    boot_est[i] <- weighted.mean(dat$x[idx], dat$w[idx])
  }
  alpha <- (1 - conf) / 2
  ci <- quantile(
    boot_est,
    probs = c(alpha, 1 - alpha),
    names = FALSE,
    na.rm = TRUE
  )
  list(
    mean  = weighted.mean(dat$x, dat$w),
    lower = ci[1],
    upper = ci[2],
    width = diff(ci),
    sig   = ifelse(ci[1] * ci[2] > 0, "*", ""),
    n_eff = n_blocks
  )
}
# -------------------------------------------------------------------------
# Spatial block bootstrap by class
# -------------------------------------------------------------------------
bootstrap_ci_by_class <- function(r_values,
                                  z_values,
                                  w_values,
                                  block_id,
                                  n_boot = 1000L,
                                  conf = 0.95) {
  # Compute spatial block bootstrap confidence intervals for weighted means of
  # r_values by class defined in z_values using weights w_values and spatial
  # blocks defined by block_id.
  if (!(
    length(r_values) == length(z_values) &&
      length(r_values) == length(w_values) &&
      length(r_values) == length(block_id)
  )) {
    stop("All inputs must have equal length.")
  }

  r_values <- as.numeric(r_values)
  w_values <- as.numeric(w_values)
  z_values <- as.character(z_values)
  block_id <- as.character(block_id)
  classes <- sort(unique(z_values[!is.na(z_values)]))
  alpha <- (1 - conf) / 2
  set.seed(26)

  results <- lapply(classes, function(cls) {
    idx <- which(
      z_values == cls &
        is.finite(r_values) &
        is.finite(w_values) &
        w_values > 0 &
        !is.na(block_id)
    )

    if (length(idx) == 0) {
      return(
        tibble::tibble(
          class = cls,
          n_pixels = 0L,
          mean_est = NA_real_,
          ci_lower = NA_real_,
          ci_upper = NA_real_,
          sig_flag = ""
        )
      )
    }

    dat <- data.frame(
      x = r_values[idx],
      w = w_values[idx],
      block = block_id[idx]
    )

    block_idx <- split(seq_len(nrow(dat)), dat$block)
    n_blocks <- length(block_idx)

    if (n_blocks < 2) {
      return(
        tibble::tibble(
          class = cls,
          n_pixels = nrow(dat),
          mean_est = weighted.mean(dat$x, dat$w),
          ci_lower = NA_real_,
          ci_upper = NA_real_,
          sig_flag = ""
        )
      )
    }

    boot_est <- numeric(n_boot)

    for (b in seq_len(n_boot)) {
      sampled_blocks <- sample(names(block_idx), n_blocks, replace = TRUE)
      idxb <- unlist(block_idx[sampled_blocks], use.names = FALSE)
      boot_est[b] <- weighted.mean(dat$x[idxb], dat$w[idxb])
    }

    ci <- quantile(boot_est,
      probs = c(alpha, 1 - alpha),
      na.rm = TRUE
    )

    tibble::tibble(
      class = as.character(cls),
      n_pixels = as.integer(nrow(dat)),
      mean_est = as.numeric(weighted.mean(dat$x, dat$w)),
      ci_lower = as.numeric(ci[1]),
      ci_upper = as.numeric(ci[2]),
      sig_flag = as.character(ifelse(ci[1] * ci[2] > 0, "*", ""))
    )
  })

  dplyr::bind_rows(results)
}
