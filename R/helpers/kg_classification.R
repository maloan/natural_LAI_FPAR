# Köppen-Geiger lookup on the existing 0.25-degree LAI grid. The nominal
# 100-arc-second classification bundled with kgc is sampled categorically at
# each LAI cell centre; climate categories are never interpolated.

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

kg_lookup <- function(xy, resolution = c("fine", "coarse")) {
  resolution <- match.arg(resolution)
  stopifnot(ncol(xy) == 2L, all(is.finite(xy)))
  if (any(xy[, 1] < -180 | xy[, 1] >= 180 |
          xy[, 2] <= -90 | xy[, 2] > 90)) {
    stop("Climate lookup coordinates are outside the global grid.")
  }

  if (resolution == "fine") {
    fine_grid <- terra::rast(
      nrows = 6480L,
      ncols = 12960L,
      xmin = -180,
      xmax = 180,
      ymin = -90,
      ymax = 90,
      crs = "EPSG:4326"
    )
    cell <- terra::cellFromXY(fine_grid, xy)
    return(as.character(kgc::getZone(kg_fine_values()[cell])))
  }

  dat <- new.env(parent = emptyenv())
  utils::data("climatezones", package = "kgc", envir = dat)
  centre <- function(x) round(x) + ifelse(x >= round(x), 0.25, -0.25)
  key <- function(lon, lat) paste(lon, lat, sep = ":")
  index <- match(
    key(centre(xy[, 1]), centre(xy[, 2])),
    key(dat$climatezones$Lon, dat$climatezones$Lat)
  )
  as.character(dat$climatezones$Cls[index])
}

kg_grid <- function(reference,
                    valid_cells,
                    resolution = "fine",
                    cache_dir) {
  stopifnot(resolution %in% c("fine", "coarse"))
  nominal_degrees <- if (resolution == "fine") 100 / 3600 else 0.5
  signature <- digest::digest(list(
    method = paste0("cell-centre-", resolution, "-v1"),
    package = as.character(utils::packageVersion("kgc")),
    dimensions = dim(reference),
    extent = as.vector(terra::ext(reference)),
    crs = terra::crs(reference),
    valid_cells = valid_cells
  ), algo = "xxhash64")
  path <- file.path(
    cache_dir,
    paste0("kg_grid_", resolution, "_", signature, ".rds")
  )
  if (file.exists(path)) return(readRDS(path))

  raw <- rep(NA_character_, terra::ncell(reference))
  raw[valid_cells] <- kg_lookup(
    terra::xyFromCell(reference, valid_cells),
    resolution = resolution
  )
  code <- raw
  code[!code %in% kg_codes()] <- NA_character_
  result <- list(
    code = code,
    raw = raw,
    signature = signature,
    nominal_degrees = nominal_degrees,
    resolution = resolution
  )
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, path)
  result
}

kg_output_dir <- function(kind = c("tables", "figures")) {
  kind <- match.arg(kind)
  here::here(
    "analysis", "results", kind,
    if (kind == "tables") "koppen_geiger" else "summaries"
  )
}
