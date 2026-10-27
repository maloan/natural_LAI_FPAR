# Shared ESA-CCI land-cover definitions for the 0.25-degree summaries.

lc_fraction_schema <- "primary22_v1"

# Preserve these 22 primary classes in the annual fraction rasters. Native
# subclasses (for example, 61 and 62) are combined with their parent class.
lc_primary_classes <- c(
  10L, 20L, 30L, 40L, 50L, 60L, 70L, 80L, 90L, 100L, 110L,
  120L, 130L, 140L, 150L, 160L, 170L, 180L, 190L, 200L, 210L, 220L
)

lc_native_to_primary <- c(
  `10` = 10L, `11` = 10L, `12` = 10L, `20` = 20L,
  `30` = 30L, `40` = 40L, `50` = 50L,
  `60` = 60L, `61` = 60L, `62` = 60L,
  `70` = 70L, `71` = 70L, `72` = 70L,
  `80` = 80L, `81` = 80L, `82` = 80L,
  `90` = 90L, `100` = 100L, `110` = 110L,
  `120` = 120L, `121` = 120L, `122` = 120L,
  `130` = 130L, `140` = 140L,
  `150` = 150L, `151` = 150L, `152` = 150L, `153` = 150L,
  `160` = 160L, `170` = 170L, `180` = 180L,
  `190` = 190L, `200` = 200L, `201` = 200L, `202` = 200L,
  `210` = 210L, `220` = 220L
)

# Fractions are summed within these 19 groups before the dominant class is
# selected. The list names are the class codes written to the dominant raster.
lc_group_members <- list(
  `10` = c(10L, 20L),
  `30` = 30L,
  `40` = 40L,
  `50` = 50L,
  `60` = 60L,
  `70` = 70L,
  `80` = 80L,
  `90` = 90L,
  `100` = c(100L, 110L),
  `120` = 120L,
  `130` = 130L,
  `140` = 140L,
  `150` = 150L,
  `160` = c(160L, 170L),
  `180` = 180L,
  `190` = 190L,
  `200` = 200L,
  `210` = 210L,
  `220` = 220L
)

lc_group_legend <- data.frame(
  lc_id = as.integer(names(lc_group_members)),
  lc_name = c(
    "Cropland",
    "Cropland-dominated mosaic",
    "Natural-vegetation-dominated mosaic",
    "Evergreen broadleaf forest",
    "Deciduous broadleaf forest",
    "Evergreen needleleaf forest",
    "Deciduous needleleaf forest",
    "Mixed forest",
    "Woody-herbaceous mosaic",
    "Shrubland",
    "Grassland",
    "Lichens and mosses",
    "Sparse vegetation",
    "Flooded tree cover",
    "Flooded shrub/herbaceous cover",
    "Urban areas",
    "Bare areas",
    "Water bodies",
    "Permanent snow and ice"
  ),
  stringsAsFactors = FALSE
)

# The paper figure includes vegetated classes and urban areas. The three
# nonvegetated groups remain in the classification and result tables.
lc_plot_ids <- setdiff(lc_group_legend$lc_id, c(200L, 210L, 220L))

lc_expected_fraction_names <- function() {
  paste0("lc_", lc_primary_classes)
}

lc_fraction_file_is_valid <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(FALSE)
  }
  x <- try(terra::rast(path), silent = TRUE)
  if (inherits(x, "try-error")) {
    return(FALSE)
  }
  identical(names(x), lc_expected_fraction_names()) &&
    terra::nrow(x) == 720L && terra::ncol(x) == 1440L &&
    terra::nlyr(x) == length(lc_primary_classes) &&
    isTRUE(all.equal(unname(as.vector(terra::ext(x))), c(-180, 180, -90, 90)))
}

lc_assert_fraction_files <- function(paths) {
  valid <- vapply(paths, lc_fraction_file_is_valid, logical(1))
  if (any(!valid)) {
    stop(
      "Invalid or missing ", lc_fraction_schema,
      " annual fraction files:\n",
      paste(paths[!valid], collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

lc_group_fraction_raster <- function(primary_fraction) {
  if (!identical(names(primary_fraction), lc_expected_fraction_names())) {
    stop("Unexpected primary land-cover fraction layers.", call. = FALSE)
  }

  grouped <- lapply(lc_group_members, function(class_ids) {
    layers <- primary_fraction[[paste0("lc_", class_ids)]]
    if (terra::nlyr(layers) == 1L) {
      return(layers)
    }
    terra::app(layers, sum, na.rm = FALSE)
  })
  grouped <- terra::rast(grouped)
  names(grouped) <- paste0("lc_", names(lc_group_members))
  grouped
}

lc_class_name <- function(class_ids) {
  lc_group_legend$lc_name[match(class_ids, lc_group_legend$lc_id)]
}
