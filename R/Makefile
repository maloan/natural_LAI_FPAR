# =============================================================================
# Makefile — natural_LAI_FPAR
# =============================================================================

R    ?= Rscript
ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)
WD   := $(ROOT)/R

RUN_TAG ?= alpha_0.1
MASKS   := CCI GLC

# The run tag selects the CCI threshold scenario.
ENV := RUN_TAG=$(RUN_TAG)

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
masks: setup
	$(call run,03_cci_mask_0p05.R)
	$(call run,05_glc_mask_0p05.R)
	$(call run,06_nonveg_static_from_cci_0p05.R)
	GRASS_SOURCE=CCI $(call run,09_luh_pasture_overlap_0p25.R)
	GRASS_SOURCE=GLC $(call run,09_luh_pasture_overlap_0p25.R)

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
# 4) Trend products (includes MK p-values via trends/compute_mk_pval.R)
# -----------------------------------------------------------------------------
.PHONY: trends trends-unmasked trends-masked
trends-unmasked:
	bash $(ROOT)/trends/01_build_georef_products.sh

trends-masked:
	@for m in $(MASKS); do \
	  echo "== TRENDS $$m (RUN_TAG=$(RUN_TAG)) ==" ; \
	  bash $(ROOT)/trends/build_trends_masked_0p25.sh $(RUN_TAG) LAI  $$m ; \
	  bash $(ROOT)/trends/build_trends_masked_0p25.sh $(RUN_TAG) FPAR $$m ; \
	done

trends: trends-unmasked trends-masked

# -----------------------------------------------------------------------------
# 5) Analyses
.PHONY: analysis
analysis: pipeline
	$(R) $(WD)/analysis/00_area_validdomain_after_nonvegetated.R
# -----------------------------------------------------------------------------
# Help
# -----------------------------------------------------------------------------
.PHONY: help
help:
	@echo "Targets:"
	@echo "  pipeline   -> setup + prep (georef+fractions+stacks) + masks + apply + agg"
	@echo "  trends     -> run unmasked + masked trend workflows (incl. MK p-values)"
	@echo "  trends-unmasked -> build georef trend products at 0.05/0.25"
	@echo "  trends-masked   -> build masked 0.25 trend products for RUN_TAG and MASKS"
	@echo "  analysis   -> pipeline + analysis/00_area_validdomain_after_nonvegetated.R"
	@echo ""
	@echo "Typical runs:"
	@echo "  make pipeline"
	@echo "  make trends"
	@echo "  make pipeline RUN_TAG=alpha_0.2  # optional sensitivity run"
