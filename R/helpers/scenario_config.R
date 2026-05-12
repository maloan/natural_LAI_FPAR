# Helper to consolidate scenario specification across analysis scripts
# Usage: scenario_spec <- create_scenario_spec()

create_scenario_spec <- function(
  cci_taus = c("tau_0.05", "tau_0.1", "tau_0.2"),
  glc_run_tag = "tau_0.1"
) {
  tibble::tibble(
    scenario = c(
      "Unmasked",
      sprintf("CCI tau=%s", sub("tau_", "", cci_taus)),
      "GLC"
    ),
    source = c("unmasked", rep("CCI", length(cci_taus)), "GLC"),
    run_tag = c(NA_character_, cci_taus, glc_run_tag)
  )
}

# Scenario order vector for factor levels
scenario_order <- function(scenario_spec) {
  as.character(scenario_spec$scenario)
}
