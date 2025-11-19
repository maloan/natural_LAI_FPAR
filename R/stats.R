## =============================================================================
# stats.R — Area-weighted statistics for binary masks (drop=1, keep=0)
#
# Purpose
#   Robust, alignment-safe evaluators for binary masks on latitude-varying
#   area grids. All ≥1 values are treated as “drop”.
#
# Functions
#   .ensure_aligned(a, b, w)  — Align b and w to a (near for masks, bilinear for w)
#   aw_mean(x, w)             — Area-weighted mean over finite support
#   p_drop(mask01, w)         — Area-weighted fraction dropped
#   pairwise_stats(a, b, w)   — Weighted contingency + Jaccard on drop class
#
# Dependencies
#   terra
## =============================================================================


# ------------------------------------------------------------------------------
# Alignment helper
# ------------------------------------------------------------------------------

.ensure_aligned <- function(a, b, w) {
  stopifnot(inherits(a, "SpatRaster"),
            inherits(b, "SpatRaster"),
            inherits(w, "SpatRaster"))

  if (!terra::compareGeom(b, a, stopOnError = FALSE)) {
    b <- terra::resample(b, a, method = "near")     # categorical
  }
  if (!terra::compareGeom(w, a, stopOnError = FALSE)) {
    w <- terra::resample(w, a, method = "bilinear") # continuous weights
  }
  list(a = a, b = b, w = w)
}


# ------------------------------------------------------------------------------
# Area-weighted mean (strict NA control)
# ------------------------------------------------------------------------------

aw_mean <- function(x, w) {
  stopifnot(inherits(x, "SpatRaster"),
            inherits(w, "SpatRaster"))

  ok  <- terra::ifel(is.finite(x) & is.finite(w), 1, NA)
  num <- terra::global(x * w * ok, "sum", na.rm = TRUE)[1, 1]
  den <- terra::global(w * ok,       "sum", na.rm = TRUE)[1, 1]

  if (is.na(den) || den == 0) NA_real_ else as.numeric(num / den)
}


# ------------------------------------------------------------------------------
# Area-weighted fraction dropped (mask ≥1 → 1)
# ------------------------------------------------------------------------------

p_drop <- function(mask01, w) {
  stopifnot(inherits(mask01, "SpatRaster"),
            inherits(w,      "SpatRaster"))

  al <- .ensure_aligned(mask01, mask01, w)
  M  <- terra::ifel(al$a >= 1, 1, 0)
  aw_mean(M, al$w)
}


# ------------------------------------------------------------------------------
# Pairwise contingency + Jaccard (drop class)
# ------------------------------------------------------------------------------

pairwise_stats <- function(a, b, w) {
  stopifnot(inherits(a, "SpatRaster"),
            inherits(b, "SpatRaster"),
            inherits(w, "SpatRaster"))

  al <- .ensure_aligned(a, b, w)
  A  <- terra::ifel(al$a >= 1, 1, 0)
  B  <- terra::ifel(al$b >= 1, 1, 0)
  W  <- al$w

  # weighted sum helper
  ws <- function(r) {
    ok <- terra::ifel(is.finite(r) & is.finite(W), 1, NA)
    terra::global(r * W * ok, "sum", na.rm = TRUE)[1, 1]
  }

  w11 <- ws((A == 1) & (B == 1))
  w10 <- ws((A == 1) & (B == 0))
  w01 <- ws((A == 0) & (B == 1))
  w00 <- ws((A == 0) & (B == 0))

  wA  <- ws(A)
  wB  <- ws(B)
  wOR <- ws((A == 1) | (B == 1))

  jacc <- if (!is.na(wOR) && wOR > 0) w11 / wOR else NA_real_

  data.frame(
    wA   = wA,
    wB   = wB,
    w00  = w00,
    w01  = w01,
    w10  = w10,
    w11  = w11,
    wAND = w11,
    wOR  = wOR,
    jacc = jacc,
    row.names = NULL
  )
}
