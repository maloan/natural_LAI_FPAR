## =============================================================================
## files.R — Filename and file-discovery helpers
## =============================================================================

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

find_one <- function(dir, pattern) {
  cand <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(cand)) {
    stop("No file matching pattern '", pattern, "' in: ", dir, call. = FALSE)
  }
  cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]
}
