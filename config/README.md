# Configuration

This folder contains one configuration file for each masking scenario.

## Configuration files

- `config_alpha_0.05.yml`
- `config_alpha_0.2.yml`
- `config_alpha_0.1.yml`
- `config_glc.yml` (GLC scenario)

The three `alpha_*` files define the CCI thresholds. `config_glc.yml` defines the GLC scenario. The `make` workflow selects each configuration in sequence.

