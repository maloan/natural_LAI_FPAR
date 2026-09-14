#!/usr/bin/env python3
"""Aggregate annual 300 m ESA-CCI land cover to 0.25-degree fractions."""

from pathlib import Path

import netCDF4
import numpy as np
import rasterio
from pyproj import Geod
from rasterio.transform import from_origin


REPO_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = REPO_DIR / "data-raw" / "ESACCI" / "ESACCI_1992-2022"
FRACTION_DIR = REPO_DIR / "analysis" / "tmp" / "lc025_fraction_yearly"
MAJORITY_DIR = REPO_DIR / "analysis" / "tmp" / "lc025_majority_yearly"

PARENT_CLASSES = np.array(
    [10, 50, 60, 70, 80, 90, 100, 120, 130, 140, 150, 160, 180, 190, 200, 210, 220],
    dtype=np.int16,
)

CLASS_TO_PARENT = {
    10: 10,
    11: 10,
    12: 10,
    20: 10,
    30: 10,
    40: 10,
    50: 50,
    60: 60,
    61: 60,
    62: 60,
    70: 70,
    71: 70,
    72: 70,
    80: 80,
    81: 80,
    82: 80,
    90: 90,
    100: 100,
    110: 100,
    120: 120,
    121: 120,
    122: 120,
    130: 130,
    140: 140,
    150: 150,
    151: 150,
    152: 150,
    153: 150,
    160: 160,
    170: 160,
    180: 180,
    190: 190,
    200: 200,
    201: 200,
    202: 200,
    210: 210,
    220: 220,
}


def input_path(year: int) -> Path:
    if year <= 2015:
        name = f"ESACCI-LC-L4-LCCS-Map-300m-P1Y-{year}-v2.0.7cds.nc"
    else:
        name = f"C3S-LC-L4-LCCS-Map-300m-P1Y-{year}-v2.1.1.nc"
    return INPUT_DIR / name


def parent_lookup() -> np.ndarray:
    class_index = {value: index + 1 for index, value in enumerate(PARENT_CLASSES)}
    lookup = np.zeros(256, dtype=np.uint8)
    for source_class, parent_class in CLASS_TO_PARENT.items():
        lookup[source_class] = class_index[parent_class]
    return lookup


def aggregate_year(year: int, lookup: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    source_path = input_path(year)
    if not source_path.exists():
        raise FileNotFoundError(f"Missing expected ESA-CCI file: {source_path}")

    fractions = np.full((len(PARENT_CLASSES), 720, 1440), np.nan, dtype=np.float32)
    target_column_key = np.repeat(
        np.arange(1440, dtype=np.int32) * len(PARENT_CLASSES),
        90,
    )
    latitude_edges = 90 - np.arange(64801, dtype=np.float64) / 360
    geod = Geod(ellps="WGS84")
    row_weight = np.empty(64800, dtype=np.float64)
    for row in range(64800):
        area, _ = geod.polygon_area_perimeter(
            [0, 1 / 360, 1 / 360, 0],
            [
                latitude_edges[row],
                latitude_edges[row],
                latitude_edges[row + 1],
                latitude_edges[row + 1],
            ],
        )
        row_weight[row] = abs(area)

    with netCDF4.Dataset(source_path) as dataset:
        source = dataset.variables["lccs_class"]
        source.set_auto_maskandscale(False)
        if source.shape != (1, 64800, 129600):
            raise ValueError(f"Unexpected lccs_class shape in {source_path}: {source.shape}")

        for source_start in range(0, 64800, 4050):
            source_stop = min(source_start + 4050, 64800)
            source_block = source[0, source_start:source_stop, :]

            for offset in range(0, source_stop - source_start, 90):
                output_row = (source_start + offset) // 90
                raw = source_block[offset : offset + 90, :]
                parent = lookup[raw]
                valid = parent > 0
                key = target_column_key[None, :] + parent.astype(np.int32) - 1
                weights = np.broadcast_to(
                    row_weight[source_start + offset : source_start + offset + 90, None],
                    raw.shape,
                )
                weighted_area = np.bincount(
                    key[valid],
                    weights=weights[valid],
                    minlength=1440 * len(PARENT_CLASSES),
                ).reshape(1440, len(PARENT_CLASSES))
                total_area = weighted_area.sum(axis=1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    fractions[:, output_row, :] = (weighted_area / total_area[:, None]).T

            print(f"  source rows {source_start + 1}-{source_stop}", flush=True)

    valid = np.isfinite(fractions).any(axis=0)
    majority_index = np.argmax(np.where(np.isfinite(fractions), fractions, -1), axis=0)
    majority = np.where(valid, PARENT_CLASSES[majority_index], -9999).astype(np.int16)
    return fractions, majority


def write_outputs(year: int, fractions: np.ndarray, majority: np.ndarray) -> None:
    transform = from_origin(-180, 90, 0.25, 0.25)
    fraction_path = FRACTION_DIR / f"lc025_fraction_{year}.tif"
    majority_path = MAJORITY_DIR / f"lc025_majority_{year}.tif"

    fraction_profile = {
        "driver": "GTiff",
        "height": 720,
        "width": 1440,
        "count": len(PARENT_CLASSES),
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": transform,
        "nodata": -9999.0,
        "compress": "deflate",
        "predictor": 3,
        "tiled": True,
        "BIGTIFF": "IF_SAFER",
    }
    with rasterio.open(fraction_path, "w", **fraction_profile) as output:
        output.write(np.where(np.isfinite(fractions), fractions, -9999.0))
        for band, class_id in enumerate(PARENT_CLASSES, start=1):
            output.set_band_description(band, f"lc_{class_id}")

    majority_profile = fraction_profile | {
        "count": 1,
        "dtype": "int16",
        "nodata": -9999,
        "predictor": 2,
    }
    with rasterio.open(majority_path, "w", **majority_profile) as output:
        output.write(majority, 1)
        output.set_band_description(1, "lc_id")


def main() -> None:
    FRACTION_DIR.mkdir(parents=True, exist_ok=True)
    MAJORITY_DIR.mkdir(parents=True, exist_ok=True)
    lookup = parent_lookup()

    for year in range(1992, 2023):
        fraction_path = FRACTION_DIR / f"lc025_fraction_{year}.tif"
        majority_path = MAJORITY_DIR / f"lc025_majority_{year}.tif"
        if fraction_path.exists() and majority_path.exists():
            print(f"Outputs exist, skipping {year}", flush=True)
            continue

        print(f"Processing {year}: {input_path(year).name}", flush=True)
        fractions, majority = aggregate_year(year, lookup)
        write_outputs(year, fractions, majority)
        print(f"Completed {year}", flush=True)


if __name__ == "__main__":
    main()
