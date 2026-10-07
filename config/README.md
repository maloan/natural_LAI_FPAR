# Configuration

This folder contains generated configuration files for the LAI and FPAR managed-land masking workflow.

Each processing scenario has its own configuration file and scripts load the file matching the active `RUN_TAG`.

## Configuration files

- `config_alpha_0.05.yml`
- `config_alpha_0.1.yml`
- `config_alpha_0.2.yml`

These files correspond to the three CCI masking thresholds used in the analysis.

All configuration files share the same structure and define:

- run metadata and time span
- input and output paths
- reference and area grids
- land-cover class mappings
- output naming conventions

## Usage

Scripts load the active configuration with:

```r
cfg <- cfg_read()
```

Run `R/00_setup.R` for each required `RUN_TAG` to generate the corresponding configuration file.

The workflow uses a 0.05° processing grid and a 0.25° analysis grid.

## Editing

Typical updates include:

- local or cluster paths
- analysis years
- land-cover class mappings
- optional region or quicklook settings