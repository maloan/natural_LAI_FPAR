# Configuration (config)

This folder holds the main settings for the natural vegetation LAI/FPAR workflow.
In practice, almost every script reads this configuration through cfg_read().

The main file is config.yml. It describes what a run should do and where it should read and write data.

## What is in this folder

- config.yml
    Active run configuration used by scripts via cfg_read().
- config_tau_0.05.yml, config_tau_0.1.yml, config_tau_0.2.yml
    Saved configuration snapshots for common tau runs.

All config files share the same schema and include:
project metadata (run tag, CRS, time span),
input/output paths,
reference and area grids,
land-cover class mappings (ESA-CCI, GLC-FCS30D, LUH2),
and output naming templates.

## What you usually edit

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

The file config.yml is generated/updated by R/00_setup.R for each run.
