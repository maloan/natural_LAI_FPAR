# =============================================================================
# Paths & tools
# =============================================================================
R    ?= Rscript
ROOT ?= $(HOME)/GitHub/natural_LAI_FPAR
WD   ?= $(ROOT)/R
RUN_TAG ?= tau_0.2
FRACTIONAL_MASKING ?= TRUE

export RUN_TAG
export FRACTIONAL_MASKING

ENV_COMMON = \
  SKIP_EXISTING=$(SKIP_EXISTING) \
  OVERWRITE=$(OVERWRITE) \
  SILENT_TIMING=$(SILENT_TIMING) \
  REMAKE_QL=$(REMAKE_QL) \
  RUN_TAG=$(RUN_TAG)

# =============================================================================
# Main target
# =============================================================================
.PHONY: all
all: full-analysis

# =============================================================================
# 0. Setup
# =============================================================================
.PHONY: config
config:
	$(ENV_COMMON) $(R) $(WD)/00_setup.R

# =============================================================================
# 1. Georeferencing & base preprocessing
# =============================================================================
.PHONY: georef
georef: config
	VAR=LAI  $(ENV_COMMON) $(R) $(WD)/01_georef_0p05.R
	VAR=FPAR $(ENV_COMMON) $(R) $(WD)/01_georef_0p05.R

.PHONY: cci-frac
cci-frac: config
	$(ENV_COMMON) $(R) $(WD)/03_cci_frac_0p05.R

.PHONY: glc-stack
glc-stack: config
	$(ENV_COMMON) $(R) $(WD)/05_glc_stack_0p05.R

# =============================================================================
# 2. Masks (CCI, GLC, LUH, abiotic)
# =============================================================================
.PHONY: cci-mask
cci-mask:
	$(ENV_COMMON) $(R) $(WD)/04_cci_mask_0p05.R

.PHONY: glc-mask
glc-mask:
	$(ENV_COMMON) $(R) $(WD)/06_glc_mask_0p05.R

.PHONY: abiotic
abiotic:
	$(ENV_COMMON) $(R) $(WD)/07_abiotic_static_from_cci_0p05.R

.PHONY: luh-overlap
luh-overlap:
	$(ENV_COMMON) $(R) $(WD)/10_luh_pasture_overlap_025.R

# =============================================================================
# 3. Apply masks (0.05°)
# =============================================================================
apply-mask-%:
	VAR=$* $(ENV_COMMON) $(R) $(WD)/11_apply_mask_0p05.R

# Convenience groups
.PHONY: apply-all
apply-all: apply-mask-LAI apply-mask-FPAR

# =============================================================================
# 4. Aggregate to 0.25°
# =============================================================================
agg-%:
	VAR=$* $(ENV_COMMON) $(R) $(WD)/12_agg_0p25.R

.PHONY: agg-all
agg-all: agg-LAI agg-FPAR

# =============================================================================
# 5. Analysis (scripts 13–17)
# =============================================================================
.PHONY: analyse-georef
analyse-georef:
	$(ENV_COMMON) $(R) $(WD)/13_analyse_georef_0p25.R

.PHONY: compare-trends
compare-trends:
	$(ENV_COMMON) $(R) $(WD)/14_compare_lai_fpar_trends.R

.PHONY: analyse-masked-trends
analyse-masked-trends:
	$(ENV_COMMON) $(R) $(WD)/15_analyse_masked_trends.R

.PHONY: compare-masked-unmasked
compare-masked-unmasked:
	$(ENV_COMMON) $(R) $(WD)/16_compare_masked_unmasked.R

.PHONY: analyse-trend-mask
analyse-trend-mask:
	$(ENV_COMMON) $(R) $(WD)/17_analyse_trend_mask.R

.PHONY: analysis
analysis: analyse-georef compare-trends analyse-masked-trends compare-masked-unmasked analyse-trend-mask

# =============================================================================
# 6. Full pipeline
# =============================================================================
.PHONY: full-pipeline
full-pipeline: georef cci-frac glc-stack cci-mask glc-mask abiotic luh-overlap apply-all agg-all

.PHONY: full-analysis
full-analysis: full-pipeline analysis

# =============================================================================
# Quicklooks
# =============================================================================
.PHONY: ql
ql:
	REMAKE_QL=TRUE $(ENV_COMMON) $(R) $(WD)/03_cci_frac_0p05.R
	REMAKE_QL=TRUE $(ENV_COMMON) $(R) $(WD)/04_cci_mask_0p05.R
	REMAKE_QL=TRUE $(ENV_COMMON) $(R) $(WD)/05_glc_stack_0p05.R
	REMAKE_QL=TRUE $(ENV_COMMON) $(R) $(WD)/06_glc_mask_0p05.R

# =============================================================================
# Help
# =============================================================================
.PHONY: help
help:
	@echo "Targets:"
	@echo "  all                 → full pipeline + all analyses"
	@echo "  georef              → georeference LAI & FPAR"
	@echo "  cci-frac, glc-stack → preprocessing"
	@echo "  cci-mask, glc-mask, abiotic, luh-overlap"
	@echo "  apply-mask-LAI / apply-mask-FPAR"
	@echo "  agg-LAI / agg-FPAR"
	@echo "  analysis            → scripts 13–17"
	@echo "  full-pipeline       → data prep + masking + aggregation"
	@echo "  full-analysis       → full pipeline + analysis"
