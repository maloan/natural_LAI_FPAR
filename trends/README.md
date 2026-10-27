# Trend Workflows

Scripts for computing annual LAI/FPAR diagnostics, OLS trends, relative trends, and Mann--Kendall significance at 0.25° resolution.

## Run

```bash
make
```

`01_build_unmasked_0p25.sh` builds the water- and ice-masked reference trends. `02_batch_build_trends_masked.sh` builds masked trends for LAI and FPAR across the three CCI thresholds and the GLC scenario. `build_trends_masked_0p25.sh` processes one scenario and variable within the batch workflow.

## Outputs

The workflows produce:

- annual mean, maximum, minimum, and amplitude
- OLS trend slopes and intercepts
- relative trends
- grid-cell Mann--Kendall p-values

Relative trends are calculated only above variable-specific mean-value thresholds to avoid unstable values near zero.

Masked outputs are stored under:

```text
output/<RUN_TAG>/eval/trend_<VAR>_<MASKTAG>/
```

Unmasked outputs are stored under:

```text
analysis/unmasked/0p25/
```

`compute_mk_pval.R` is called internally by the trend workflows to calculate grid-cell Mann--Kendall p-values.
