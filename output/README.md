# **Final Outputs**

This directory holds all generated products from the LAI/FPAR natural-vegetation workflow.
Outputs are organized by *run tag* (e.g., `tau_0.05`, `tau_0.1`, `tau_0.2`), each reflecting a particular masking configuration.

## **Structure**

Each run tag contains:

### **1. `masked_0p05/`**

Masked monthly LAI/FPAR fields at native **0.05°** resolution.

* Values reflect the selected mask (CCI, GLC, LUH), applied per month.
* Includes optional abiotic and LUH-overlap overlays depending on environment flags.
* Quicklooks (PNG) are generated for selected diagnostic months.

### **2. `masked_0p25/`**

Area-weighted **0.05° → 0.25°** aggregates.

* Monthly 0.25° LAI/FPAR rasters.
* Retain provenance tags and GDAL metadata.
* Used by the trend and AOI analysis scripts in `analysis/`.

### **3. `masks/`**

Derived binary drop-masks (1 = drop, 0 = keep) for:

* **CCI used-mask** (τ, k thresholds, band selection)
* **GLC used ≥ N years**
* **LUH pasture/grass overlap** (Gmin/Pmin/α rules)
* **Static abiotic** masks (water/ice/bare)

Masks are stored at both **0.05°** and replicated **0.25°**.

## **Notes**

* All raster outputs are float or byte **GeoTIFFs** with standardized GDAL compression options.
* Each run tag includes quicklook previews and optional manifest CSV summaries.
