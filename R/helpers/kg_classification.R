# Köppen-Geiger lookup on the existing 0.25-degree LAI grid. Climate
# categories are assigned categorically and are never interpolated.

kg_codes <- function() kgc::getZone(seq_len(31L))

kg_fine_values <- local({
  values <- NULL
  function() {
    if (is.null(values)) {
      dat <- new.env(parent = emptyenv())
      utils::data("kmz", package = "kgc", envir = dat)
      values <<- dat$kmz
    }
    values
  }
})

kg_dominant_index <- function(area_by_class, centre_class) {
  stopifnot(
    is.matrix(area_by_class),
    length(centre_class) == ncol(area_by_class)
  )

  maximum <- apply(area_by_class, 2L, max)
  dominant <- max.col(t(area_by_class), ties.method = "first")

  # Exact boundary ties are resolved using the central source cell. This
  # avoids assigning a class from the arbitrary ordering of class codes when
  # the centre class is one of the equally extensive classes.
  valid_centre <- is.finite(centre_class) &
    centre_class >= 1L & centre_class <= nrow(area_by_class)
  centre_area <- rep(NA_real_, ncol(area_by_class))
  centre_area[valid_centre] <- area_by_class[cbind(
    centre_class[valid_centre],
    which(valid_centre)
  )]
  tolerance <- sqrt(.Machine$double.eps) * pmax(1, maximum)
  centre_is_tied <- valid_centre & abs(centre_area - maximum) <= tolerance
  dominant[centre_is_tied] <- centre_class[centre_is_tied]
  dominant[!is.finite(maximum) | maximum <= 0] <- NA_integer_
  dominant
}

kg_grid_dominant <- function(reference, valid_cells, cache_dir) {
  expected <- terra::rast(
    nrows = 720L,
    ncols = 1440L,
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 90,
    crs = "EPSG:4326"
  )
  terra::compareGeom(reference, expected, stopOnError = TRUE)

  valid_cells <- sort(unique(as.integer(valid_cells)))
  if (!length(valid_cells) ||
      anyNA(valid_cells) ||
      any(valid_cells < 1L | valid_cells > terra::ncell(reference))) {
    stop("valid_cells must contain valid 0.25-degree grid-cell indices.")
  }

  signature <- digest::digest(list(
    method = "area-dominant-9x9-centre-tie-v1",
    kgc_version = as.character(utils::packageVersion("kgc")),
    dimensions = dim(reference),
    extent = as.vector(terra::ext(reference)),
    crs = terra::crs(reference),
    valid_cells = valid_cells
  ), algo = "xxhash64")

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(
    cache_dir,
    paste0("kg_grid_dominant_", signature, ".rds")
  )
  if (file.exists(path)) return(readRDS(path))

  fine_values <- as.integer(kg_fine_values())
  codes3 <- kg_codes()
  groups3 <- substr(codes3, 1L, 2L)
  codes2 <- unique(groups3)
  group_index <- match(groups3, codes2)
  n_class3 <- length(codes3)
  n_col <- terra::ncol(reference)
  source_n_col <- 9L * n_col
  source_cols <- rep(seq_len(n_col), each = 9L)

  fine_strip <- terra::rast(
    nrows = 6480L,
    ncols = 1L,
    xmin = 0,
    xmax = 1 / 36,
    ymin = -90,
    ymax = 90,
    crs = "EPSG:4326"
  )
  source_row_area <- terra::values(
    terra::cellSize(fine_strip, unit = "km"),
    mat = FALSE
  )

  result3 <- rep(NA_character_, terra::ncell(reference))
  result2 <- rep(NA_character_, terra::ncell(reference))

  for (target_row in seq_len(terra::nrow(reference))) {
    area3 <- matrix(0, nrow = n_class3, ncol = n_col)
    centre3 <- rep(NA_integer_, n_col)

    for (row_offset in seq_len(9L)) {
      source_row <- (target_row - 1L) * 9L + row_offset
      first <- (source_row - 1L) * source_n_col
      row_values <- fine_values[first + seq_len(source_n_col)]
      valid <- !is.na(row_values) &
        row_values >= 1L & row_values <= n_class3

      bins <- row_values[valid] +
        (source_cols[valid] - 1L) * n_class3
      area3 <- area3 + matrix(
        tabulate(bins, nbins = n_class3 * n_col),
        nrow = n_class3
      ) * source_row_area[source_row]

      if (row_offset == 5L) {
        centre3 <- row_values[seq(5L, source_n_col, by = 9L)]
        centre3[centre3 < 1L | centre3 > n_class3] <- NA_integer_
      }
    }

    area2 <- rowsum(area3, group_index, reorder = FALSE)
    centre2 <- group_index[centre3]
    dominant3 <- kg_dominant_index(area3, centre3)
    dominant2 <- kg_dominant_index(area2, centre2)

    cells <- (target_row - 1L) * n_col + seq_len(n_col)
    result3[cells] <- codes3[dominant3]
    result2[cells] <- codes2[dominant2]
  }

  outside <- rep(TRUE, terra::ncell(reference))
  outside[valid_cells] <- FALSE
  result3[outside] <- NA_character_
  result2[outside] <- NA_character_

  result <- list(
    code = result3,
    code2 = result2,
    signature = signature,
    nominal_degrees = 100 / 3600,
    resolution = "area-dominant-9x9"
  )

  temporary <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = cache_dir
  )
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(result, temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not move the completed Köppen-Geiger cache into place.")
  }
  result
}

kg_output_dir <- function(kind = c("tables", "figures")) {
  kind <- match.arg(kind)
  here::here(
    "analysis", "results", kind,
    if (kind == "tables") "koppen_geiger" else "summaries"
  )
}
