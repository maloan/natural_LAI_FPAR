# =============================================================================
# Makefile — natural_LAI_FPAR
# =============================================================================

R       ?= Rscript
ROOT    ?= $(HOME)/GitHub/natural_LAI_FPAR
WD      := $(ROOT)/R

RUN_TAG ?= tau_0.1
TAU_CCI ?= 0.1
MASKS   ?= CCI GLC
SKIP_EXISTING?= TRUE
REMAKE_QL?= FALSE

# Minimal environment passed to all scripts
ENV := RUN_TAG=$(RUN_TAG) TAU_CCI=$(TAU_CCI) REMAKE_QL=$(REMAKE_QL) SKIP_EXISTING=$(SKIP_EXISTING)

# Convenience macro
run = $(ENV) $(R) $(WD)/$1

.PHONY: all
all: analysis

# -----------------------------------------------------------------------------
# 0) Setup
# -----------------------------------------------------------------------------
.PHONY: setup
setup:
	$(call run,00_setup.R)

# -----------------------------------------------------------------------------
# 1) Georef + fractions / stacks
# -----------------------------------------------------------------------------
.PHONY: georef cci-frac glc-stack prep
georef: setup
	VAR=LAI  $(call run,01_georef_0p05.R)
	VAR=FPAR $(call run,01_georef_0p05.R)

cci-frac: setup
	$(call run,02_cci_frac_0p05.R)

glc-stack: setup
	$(call run,04_glc_stack_0p05.R)

prep: georef cci-frac glc-stack

# -----------------------------------------------------------------------------
# 2) Masks (build all base masks)
# -----------------------------------------------------------------------------
.PHONY: masks
masks:
	$(call run,03_cci_mask_0p05.R)
	$(call run,05_glc_mask_0p05.R)
	$(call run,06_nonveg_static_from_cci_0p05.R)
	$(call run,09_luh_pasture_overlap_0p25.R)

# -----------------------------------------------------------------------------
# 3) Apply masks (0.05°) + aggregate (0.25°)
# -----------------------------------------------------------------------------
.PHONY: apply agg
apply:
	@for m in $(MASKS); do \
	  echo "== APPLY $$m ==" ; \
	  MASK=$$m VAR=LAI  $(call run,10_apply_mask_0p05.R) ; \
	  MASK=$$m VAR=FPAR $(call run,10_apply_mask_0p05.R) ; \
	done

agg:
	@for m in $(MASKS); do \
	  echo "== AGG $$m ==" ; \
	  MASK=$$m VAR=LAI  $(call run,11_agg_0p5.R) ; \
	  MASK=$$m VAR=FPAR $(call run,11_agg_0p5.R) ; \
	  MASK=$$m VAR=LAI  $(call run,11_agg_0p25.R) ; \
	  MASK=$$m VAR=FPAR $(call run,11_agg_0p25.R) ; \
	done

# Full data pipeline
.PHONY: pipeline
pipeline: prep masks apply agg

# -----------------------------------------------------------------------------
# 4) Analyses
analysis: pipeline
.PHONY: analysis
analysis: pipeline
	$(R) $(WD)/analysis/00_area_validdomain_after_nonvegetated.R
# -----------------------------------------------------------------------------
# Quicklooks (explicit; no extra switches)
# -----------------------------------------------------------------------------
.PHONY: ql
ql:
	REMAKE_QL=TRUE $(call run,02_cci_frac_0p05.R)
	REMAKE_QL=TRUE $(call run,03_cci_mask_0p05.R)
	REMAKE_QL=TRUE $(call run,04_glc_stack_0p05.R)
	REMAKE_QL=TRUE $(call run,05_glc_mask_0p05.R)

# -----------------------------------------------------------------------------
# Help
# -----------------------------------------------------------------------------
.PHONY: help
help:
	@echo "Targets:"
	@echo "  pipeline   -> setup + prep (georef+fractions+stacks) + masks + apply + agg"
	@echo "  analysis   -> pipeline (analysis scripts are in R/analysis folder)"
	@echo "  ql         -> regenerate quicklooks for key steps"
	@echo ""
	@echo "Typical runs:"
	@echo "  make pipeline RUN_TAG=tau_0.1"
	@echo "  make pipeline MASKS='CCI'"
	@echo "  make pipeline MASKS='GLC'"
