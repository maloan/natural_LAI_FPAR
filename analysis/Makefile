# ==== paths & tools ====
R              ?= Rscript
ROOT           ?= $(HOME)/GitHub/natural_LAI_FPAR
WD             ?= $(ROOT)/src
RUN_TAG        ?= baseline
COMMON_OUT     ?= $(ROOT)/data/common

# ==== environment flags ====
SKIP_EXISTING      ?= FALSE
OVERWRITE          ?= TRUE
SILENT_TIMING      ?= TRUE
REMAKE_QL          ?= TRUE

# ==== mask settings ====
USED_N_YEARS       ?= 3
TAU_CCI            ?= 0.1
K_CCI              ?= 3
CCI_BAND           ?= frac_fused

TAU_WATER          ?= 0.05
TAU_ICE            ?= 0.05
TAU_BARE           ?= 0.30
ABIOTIC_Y1         ?= 1992
ABIOTIC_Y2         ?= 2020
APPLY_ABIOTIC_STATIC ?= TRUE

GRASS_SOURCE       ?= CCI
G_MIN              ?= 0.05
P_MIN              ?= 0.20
ALPHA              ?= 0.50
LUH_AVG_START      ?= 1992
LUH_AVG_END        ?= 2020
APPLY_LUH_OVERLAP  ?= TRUE

ENV_COMMON = SKIP_EXISTING=$(SKIP_EXISTING) OVERWRITE=$(OVERWRITE) SILENT_TIMING=$(SILENT_TIMING) REMAKE_QL=$(REMAKE_QL)

# ==== master target ====
.PHONY: all
all: all-full

.PHONY: all-full
all-full: precompute \
          apply-lai-cci-luhov apply-fpar-cci-luhov \
          apply-lai-glc-luhov apply-fpar-glc-luhov \
          agg-lai-cci-luhov agg-fpar-cci-luhov \
          agg-lai-glc-luhov agg-fpar-glc-luhov \
          eval-masks eval-outputs

# ==== preprocessing ====
.PHONY: precompute
precompute: config georef-lai georef-fpar cci-frac glc-yearstack

.PHONY: config
config:
	$(R) $(WD)/00_setup.R

.PHONY: georef-lai
georef-lai:
	VAR=LAI $(ENV_COMMON) $(R) $(WD)/01_georef_0p05.R

.PHONY: georef-fpar
georef-fpar:
	VAR=FPAR $(ENV_COMMON) $(R) $(WD)/01_georef_0p05.R

.PHONY: cci-frac
cci-frac:
	$(ENV_COMMON) $(R) $(WD)/02_cci_frac_0p05.R

.PHONY: glc-yearstack
glc-yearstack:
	$(ENV_COMMON) $(R) $(WD)/03_glc_stack_0p05.R

# ==== masks ====
.PHONY: cci-mask
cci-mask: cci-frac
	CCI_BAND=$(CCI_BAND) TAU_CCI=$(TAU_CCI) K_CCI=$(K_CCI) \
	$(ENV_COMMON) $(R) $(WD)/04_cci_mask_0p05.R

.PHONY: glc-mask
glc-mask: glc-yearstack
	USED_N_YEARS=$(USED_N_YEARS) $(ENV_COMMON) $(R) $(WD)/04_glc_mask_0p05.R

.PHONY: luh-overlap-mask
luh-overlap-mask:
	GRASS_SOURCE=$(GRASS_SOURCE) G_MIN=$(G_MIN) P_MIN=$(P_MIN) ALPHA=$(ALPHA) \
	LUH_AVG_START=$(LUH_AVG_START) LUH_AVG_END=$(LUH_AVG_END) \
	$(ENV_COMMON) $(R) $(WD)/04_luh_pasture_overlap_025.R

.PHONY: abiotic-cci
abiotic-cci:
	TAU_WATER=$(TAU_WATER) TAU_ICE=$(TAU_ICE) TAU_BARE=$(TAU_BARE) \
	ABIOTIC_Y1=$(ABIOTIC_Y1) ABIOTIC_Y2=$(ABIOTIC_Y2) \
	$(ENV_COMMON) $(R) $(WD)/05_abiotic_static_from_cci_0p05.R

# ==== apply masks (with LUH overlay) ====
define apply_template
.PHONY: apply-$(1)-$(2)-luhov
apply-$(1)-$(2)-luhov: georef-$(1) $(2)-mask luh-overlap-mask abiotic-cci
	VAR=$(1) MASK=$(2) $(if $(filter CCI,$(2)),CCI_BAND=$(CCI_BAND)) \
	USED_N_YEARS=$(USED_N_YEARS) \
	APPLY_ABIOTIC_STATIC=$(APPLY_ABIOTIC_STATIC) \
	APPLY_LUH_OVERLAP=TRUE \
	GRASS_SOURCE=$(GRASS_SOURCE) G_MIN=$(G_MIN) P_MIN=$(P_MIN) ALPHA=$(ALPHA) \
	LUH_AVG_START=$(LUH_AVG_START) LUH_AVG_END=$(LUH_AVG_END) \
	$(ENV_COMMON) $(R) $(WD)/07_apply_mask_0p05.R
endef

$(foreach var,LAI FPAR, $(foreach m,CCI GLC, $(eval $(call apply_template,$(var),$(m)))))

# ==== aggregation ====
define agg_template
.PHONY: agg-$(1)-$(2)-luhov
agg-$(1)-$(2)-luhov: apply-$(1)-$(2)-luhov
	VAR=$(1) MASK=$(2) $(ENV_COMMON) $(R) $(WD)/08_agg_0p25.R
endef

$(foreach var,LAI FPAR, $(foreach m,CCI GLC, $(eval $(call agg_template,$(var),$(m)))))

# ==== evaluation ====
.PHONY: eval-masks
eval-masks:
	$(ENV_COMMON) $(R) $(WD)/10_eval_masks.R

.PHONY: eval-outputs
eval-outputs:
	$(R) $(WD)/11_eval_outputs.R
