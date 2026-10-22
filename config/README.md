# Configuration

This folder contains generated configuration files for the LAI and FPAR managed-land masking workflow.

Each processing scenario has its own configuration file and scripts load the file matching the active `RUN_TAG`.

## Configuration files

- `config_alpha_0.05.yml`
- `config_alpha_0.2.yml`
- `config_alpha_0.1.yml`
- `config_glc.yml` (GLC scenario, with no CCI managed-land alpha)

The three `alpha_*` files correspond to CCI masking thresholds. `config_glc.yml`
selects the separate GLC scenario.

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
To regenerate configuration without rewriting reference grids, use
`CONFIG_ONLY=1 RUN_TAG=glc Rscript R/00_setup.R` (or a CCI alpha tag).

The workflow uses a 0.05° processing grid and a 0.25° analysis grid.

