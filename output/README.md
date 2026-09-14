# Final Outputs (output)

This folder contains all generated products from the natural vegetation LAI/FPAR workflow. Outputs are organized by run tag (`alpha_0.05`, `alpha_0.1`, and `alpha_0.2`), where each tag represents a CCI masking threshold. Downstream analysis and figures are built from files in this directory.

## Folder layout

```text
output/
├── alpha_0.05/
├── alpha_0.1/
└── alpha_0.2/
```

Run-tag folders are the canonical output namespaces used by downstream analyses. Unmasked baseline georeferenced products are stored under `analysis/unmasked/`. The GLC branch is stored only under `alpha_0.1` but its mask is independent of the CCI threshold.

## Files created for later research:

`output/chapter2/` contains:

- `fpar_unmasked_0p5_monthly_1982-2024.nc`: unmasked monthly fAPAR, aggregated from 0.05 to 0.5 degree using grid-cell area weights.
- `mask_CCI_alpha0p10_pasture_any_0p5.nc`: unified binary exclusion mask using only CCI alpha 0.1 and pasture (1=drop, 0=keep).
- `mask_CCI_alpha0p10_pasture_excluded_fraction_0p5.tif`: fraction of each 0.5-degree cell excluded by the two masks.

## Run-tag folder content

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
- mask_luh_overlap: Pasture-grass overlap masks.
- mask_nonvegetated: Static water and ice masks.

Most masks are generated at 0.05 degree and propagated to 0.25 degree when required.

### eval

Evaluation and trend products derived from masked 0.25 degree fields.

- Trend maps and summary statistics for LAI and FPAR.
- Separate outputs by masking source (CCI and GLC).
- Direct inputs to manuscript figures and tables.
