# Final Outputs (output)

This folder contains all generated products from the natural vegetation LAI/FPAR workflow.
Outputs are organized by run tag (for example tau_0.05, tau_0.1, tau_0.2), where each run tag represents a specific masking threshold.

Downstream analysis and figures are built from files in this directory.

## Top-Level Layout

```text
output/
├── tau_0.05/
├── tau_0.1/
└── tau_0.2/
```

Run-tag folders are the canonical output namespaces used by downstream analyses.
Unmasked baseline georeferenced products are stored under `analysis/unmasked/`.

## What Is Inside Each tau_<run_tag>

Each run-tag directory contains the core processing outputs:

### masked_0p05

Masked monthly LAI and FPAR at native 0.05 degree resolution.

- Separate LAI and FPAR subfolders.
- Masks reflect the selected setup (CCI or GLC), static non-vegetated masks and LUH pasture/grass overlap filtering mask.

These files are the immediate inputs to aggregation.

### masked_0p25

Area-weighted monthly aggregates on the 0.25 degree grid.

- Aggregated from 0.05 degree products.
- Used for trend, zonal, and global analyses.

### masks

Binary masks used in the workflow (1 = drop, 0 = keep):

- mask_cci: CCI fractional used-land masks.
- mask_glc: GLC majority used-land masks.
- mask_luh: LUH pasture/rangeland masks.
- mask_luh_overlap: Pasture-grass overlap masks.
- mask_nonvegetated: Static water and ice masks.

Most masks are generated at 0.05 degree and propagated to 0.25 degree when required.

### eval

Evaluation and trend products derived from masked 0.25 degree fields.

- Trend maps and summary statistics for LAI and FPAR.
- Separate outputs by masking source (CCI and GLC).
- Direct inputs to manuscript figures and tables.
