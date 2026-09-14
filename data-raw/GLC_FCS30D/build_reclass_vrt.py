#!/usr/bin/env python3
"""Build a temporary VRT that reclassifies GLC pixels while they are read."""

import argparse
import xml.etree.ElementTree as ET

from osgeo import gdal

gdal.UseExceptions()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("output")
    parser.add_argument("kind", choices=("class", "grass"))
    args = parser.parse_args()

    dataset = gdal.Open(args.source)
    if dataset is None:
        raise RuntimeError(f"Could not open {args.source}")

    translated = gdal.Translate(args.output, dataset, format="VRT")
    if translated is None:
        raise RuntimeError(f"Could not create {args.output}")
    translated = None
    dataset = None

    if args.kind == "class":
        expression = "np.where((a == 0) | (a == 250), 255, a)"
    else:
        expression = "np.where((a == 0) | (a == 250), 255, a == 130)"

    code = f"""
import numpy as np

def reclass(in_ar, out_ar, *args, **kwargs):
    a = in_ar[0]
    out_ar[:] = {expression}
"""

    tree = ET.parse(args.output)
    root = tree.getroot()
    for band in root.findall("VRTRasterBand"):
        band.set("subClass", "VRTDerivedRasterBand")
        band.insert(0, ET.Element("NoDataValue"))
        band[0].text = "255"
        band.insert(1, ET.Element("PixelFunctionType"))
        band[1].text = "reclass"
        band.insert(2, ET.Element("PixelFunctionLanguage"))
        band[2].text = "Python"
        band.insert(3, ET.Element("PixelFunctionCode"))
        band[3].text = code

    tree.write(args.output, encoding="UTF-8", xml_declaration=True)


if __name__ == "__main__":
    main()
