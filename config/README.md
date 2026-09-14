# Configuration (config)

This folder holds the main settings for the natural vegetation LAI/FPAR workflow. Almost every script reads this configuration through cfg_read(). Each scenario has its own generated configuration file. Scripts require the exact file matching `RUN_TAG`.

## Folder content

- config_alpha_0.05.yml, config_alpha_0.1.yml, config_alpha_0.2.yml
-> Generated configurations for the three CCI masking thresholds.

All config files share the same schema and include:
project metadata (run tag, CRS, time span), input/output paths, reference and area grids, land-cover class mappings (ESA-CCI, GLC-FCS30D, LUH2), and output naming templates.

## To edit

Most updates are small and focused:

- Paths to local or cluster data locations.
- Year windows for LAI/FPAR or land-cover inputs.
- Class mappings if source products change.
- Region or quicklook settings, if needed.

## How it is used

Scripts load the configuration like this:

```r
cfg <- cfg_read()
```

Run `R/00_setup.R` once for each required `RUN_TAG` to create the corresponding configuration file. Only the 0.05-degree processing grid and 0.25-degree analysis grid are created.
