# Trend Workflows

Scripts for computing annual LAI/FPAR diagnostics, OLS trends, relative trends, and Mann--Kendall significance at 0.25° resolution.

## Main scripts

**Unmasked baseline**

```bash
./01_build_unmasked_0p25.sh
```

Builds the baseline used as `unmasked` in the analysis. Only the common water and permanent snow or ice mask is applied.

**Masked trends**

```bash
./build_trends_masked_0p25.sh alpha_0.1 LAI CCI
```

Arguments:

- `ALPHA`: run directory, for example `alpha_0.1`
- `VAR`: `LAI` or `FPAR`
- `MASKTAG`: `CCI` or `GLC`

**Batch processing**

```bash
./02_batch_build_trends_masked.sh
```

Runs the three CCI thresholds `0.05`, `0.10`, and `0.20` for LAI and FPAR. The GLC branch is run once for each variable and stored under `alpha_0.1`. The GLC mask itself does not depend on the CCI threshold.

## Outputs

The workflows produce:

- annual mean, maximum, minimum, and amplitude
- OLS trend slopes and intercepts
- relative trends
- grid-cell Mann--Kendall p-values

Relative trends are calculated only above variable-specific mean-value thresholds to avoid unstable values near zero.

Masked outputs are stored under:

```text
output/<ALPHA>/eval/trend_<VAR>_<MASKTAG>/
```

Unmasked outputs are stored under:

```text
analysis/unmasked/0p25/
```

`compute_mk_pval.R` is called internally by the trend workflows to calculate grid-cell Mann--Kendall p-values.