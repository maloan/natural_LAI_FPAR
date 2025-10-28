## =============================================================================
# stats.R — Area-weighted stats for binary masks (drop=1, keep=0)
#
# Purpose
#   Provide robust, alignment-safe evaluators for comparing and summarizing
#   binary masks (1=drop, 0=keep) on a latitude-varying area grid.
#
# Functions
#   .ensure_aligned(a, b, w) — Align B and weights W to A’s grid (near/bilinear)
#   aw_mean(x, w)            — Area-weighted mean over finite support
#   p_drop(mask01, w)        — Area-weighted fraction dropped (mask01 ≥1 → 1)
#   pairwise_stats(a, b, w)  — Area-weighted contingency + Jaccard on drop=1
#
# Inputs
#   - a, b     : terra::SpatRaster binary masks (1=drop, 0=keep; any ≥1 treated as 1)
#   - mask01   : same as above, single mask
#   - w        : terra::SpatRaster of cell areas (km²) on some lon/lat grid
#
# Outputs
#   - p_drop() : numeric scalar, fraction of *area* flagged as drop
#   - pairwise_stats(): data.frame with weighted w00, w01, w10, w11, wAND, wOR,
#                       marginals (wA, wB), and Jaccard index on drop class
#
# Dependencies
#   Packages: terra
#
# Environment
#   - Optionally expects a global `area005` weights raster if `w` is not provided
#     by the caller (callers should pass `w` explicitly for clarity).
#
# Notes
#   - Resampling policy: nearest for categorical masks; bilinear for weights.
#   - Finite-support masking ensures NA/Inf in any input does not bias results.
#   - All logical thresholds are explicit (mask >= 1 → 1; else 0).
## =============================================================================

# --- Alignment helper (local, zero-dependency on other helpers) ---------------
.ensure_aligned <- function(a, b, w) {
  stopifnot(inherits(a, "SpatRaster"),
            inherits(b, "SpatRaster"),
            inherits(w, "SpatRaster"))
  if (!terra::compareGeom(b, a, stopOnError = FALSE)) {
    b <- terra::resample(b, a, method = "near")       # categorical
  }
  if (!terra::compareGeom(w, a, stopOnError = FALSE)) {
    w <- terra::resample(w, a, method = "bilinear")   # continuous weights
  }
  list(a = a, b = b, w = w)
}

# --- Area-weighted mean with strict NA control --------------------------------
aw_mean <- function(x, w) {
  stopifnot(inherits(x, "SpatRaster"), inherits(w, "SpatRaster"))
  ok  <- terra::ifel(is.finite(x) & is.finite(w), 1, NA)
  num <- terra::global(x * w * ok, "sum", na.rm = TRUE)[1, 1]
  den <- terra::global(w * ok, "sum", na.rm = TRUE)[1, 1]
  if (is.na(den) || den == 0) {
    return(NA_real_)
  } else {
    return(as.numeric(num / den))
  }
}

# --- Fraction of dropped area (1=drop, 0=keep) --------------------------------
p_drop <- function(mask01, w) {
  stopifnot(inherits(mask01, "SpatRaster"), inherits(w, "SpatRaster"))
  al <- .ensure_aligned(mask01, mask01, w)  # align weights to mask grid
  M  <- terra::ifel(al$a >= 1, 1, 0)
  return(aw_mean(M, al$w))
}

# --- Pairwise contingency + Jaccard (on "drop" class) -------------------------
pairwise_stats <- function(a, b, w) {
  stopifnot(inherits(a, "SpatRaster"),
            inherits(b, "SpatRaster"),
            inherits(w, "SpatRaster"))

  al <- .ensure_aligned(a, b, w)
  A  <- terra::ifel(al$a >= 1, 1, 0)
  B  <- terra::ifel(al$b >= 1, 1, 0)
  W  <- al$w

  ws <- function(r) {
    ok <- terra::ifel(is.finite(r) & is.finite(W), 1, NA)
    terra::global(r * W * ok, "sum", na.rm = TRUE)[1, 1]
  }

  w11 <- ws(terra::ifel((A == 1) & (B == 1), 1, 0))
  w10 <- ws(terra::ifel((A == 1) & (B == 0), 1, 0))
  w01 <- ws(terra::ifel((A == 0) & (B == 1), 1, 0))
  w00 <- ws(terra::ifel((A == 0) & (B == 0), 1, 0))
  wA  <- ws(A)
  wB  <- ws(B)
  wOR <- ws(terra::ifel((A == 1) | (B == 1), 1, 0))
  jacc <- if (!is.na(wOR) &&
              wOR > 0) {
    w11 / wOR
  } else {
    NA_real_
  }

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
