# ==== paths & tools ====
R ?= Rscript
ROOT ?= $(HOME)/GitHub/natural_LAI_FPAR
WD ?= $(ROOT)/R
RUN_TAG ?= tau_0.2
export RUN_TAG

# ==== env toggles ====
SKIP_EXISTING ?= FALSE
OVERWRITE     ?= TRUE
SILENT_TIMING ?= FALSE
REMAKE_QL     ?= TRUE

# GLC mask
USED_N_YEARS  ?= 3

# CCI mask
TAU_CCI  ?= 0.2
K_CCI    ?= 3
CCI_BAND ?= frac_fused

# Abiotic (CCI)
TAU_WATER ?= 0.05
TAU_ICE   ?= 0.05
TAU_BARE  ?= 0.30
ABIOTIC_Y1   ?= 1992
ABIOTIC_Y2   ?= 2020
APPLY_ABIOTIC_STATIC ?= TRUE

# LUH overlap (0.25°, replicate to 0.05°)
GRASS_SOURCE  ?= CCI
G_MIN         ?= 0.1
P_MIN         ?= 0.1
ALPHA         ?= 0.50
LUH_AVG_START ?= 1992
LUH_AVG_END   ?= 2020
APPLY_LUH_OVERLAP ?= TRUE

# Common env
ENV_COMMON = SKIP_EXISTING=$(SKIP_EXISTING) OVERWRITE=$(OVERWRITE) SILENT_TIMING=$(SILENT_TIMING) REMAKE_QL=$(REMAKE_QL) RUN_TAG=$(RUN_TAG)

# ==== default ====
.PHONY: all
all: lai-cci fpar-cci

# ---- optional meta-targets ----
.PHONY: all-cci-glc
all-cci-glc: lai-cci fpar-cci lai-glc fpar-glc

.PHONY: all-full
all-full: luh-overlap-mask lai-cci-luhov fpar-cci-luhov lai-glc-luhov fpar-glc-luhov



# ==== precomputable shared steps ====
.PHONY: precompute
precompute: georef-lai georef-fpar cci-frac glc-yearstack

.PHONY: precompute-quick
precompute-quick: config cci-frac glc-yearstack

# ==== core steps ====
.PHONY: config
config:
	$(ENV_COMMON) $(R) $(WD)/00_setup.R

.PHONY: georef-lai
georef-lai: config
	VAR=LAI  $(ENV_COMMON) $(R) $(WD)/01_georef_0p05.R

.PHONY: georef-fpar
georef-fpar: config
	VAR=FPAR $(ENV_COMMON) $(R) $(WD)/01_georef_0p05.R

.PHONY: cci-frac
cci-frac: config
	$(ENV_COMMON) $(R) $(WD)/02_cci_frac_0p05.R

.PHONY: abiotic-cci
abiotic-cci: config
	TAU_WATER=$(TAU_WATER) TAU_ICE=$(TAU_ICE) TAU_BARE=$(TAU_BARE) \
	ABIOTIC_Y1=$(ABIOTIC_Y1) ABIOTIC_Y2=$(ABIOTIC_Y2) \
	$(ENV_COMMON) $(R) $(WD)/05_abiotic_static_from_cci_0p05.R

.PHONY: cci-mask
cci-mask:
	CCI_BAND=$(CCI_BAND) TAU_CCI=$(TAU_CCI) K_CCI=$(K_CCI) \
	$(ENV_COMMON) $(R) $(WD)/04_cci_mask_0p05.R

.PHONY: glc-yearstack
glc-yearstack: config
	$(ENV_COMMON) $(R) $(WD)/03_glc_stack_0p05.R

.PHONY: glc-mask
glc-mask:
	USED_N_YEARS=$(USED_N_YEARS) $(ENV_COMMON) $(R) $(WD)/04_glc_mask_0p05.R

.PHONY: luh-use
luh-use: config
	$(ENV_COMMON) $(R) $(WD)/07_luh_use_masks.R

.PHONY: luh-overlap-mask
luh-overlap-mask: config cci-frac
	GRASS_SOURCE=$(GRASS_SOURCE) G_MIN=$(G_MIN) P_MIN=$(P_MIN) ALPHA=$(ALPHA) \
	LUH_AVG_START=$(LUH_AVG_START) LUH_AVG_END=$(LUH_AVG_END) \
	$(ENV_COMMON) $(R) $(WD)/04_luh_pasture_overlap_025.R

# ==== apply masks @0.05° ====
.PHONY: apply-lai-cci
apply-lai-cci: cci-mask abiotic-cci
	VAR=LAI  MASK=CCI  CCI_BAND=$(CCI_BAND) TAU_CCI=$(TAU_CCI) \
	APPLY_ABIOTIC_STATIC=$(APPLY_ABIOTIC_STATIC) \
	$(ENV_COMMON) $(R) $(WD)/07_apply_mask_0p05.R

.PHONY: apply-fpar-cci
apply-fpar-cci: cci-mask abiotic-cci
	VAR=FPAR MASK=CCI  CCI_BAND=$(CCI_BAND) TAU_CCI=$(TAU_CCI) \
	APPLY_ABIOTIC_STATIC=$(APPLY_ABIOTIC_STATIC) \
	$(ENV_COMMON) $(R) $(WD)/07_apply_mask_0p05.R

.PHONY: apply-lai-glc
apply-lai-glc: glc-mask abiotic-cci
	VAR=LAI  MASK=GLC  USED_N_YEARS=$(USED_N_YEARS) TAU_CCI=$(TAU_CCI) \
	APPLY_ABIOTIC_STATIC=$(APPLY_ABIOTIC_STATIC) \
	$(ENV_COMMON) $(R) $(WD)/07_apply_mask_0p05.R

.PHONY: apply-fpar-glc
apply-fpar-glc: glc-mask abiotic-cci
	VAR=FPAR MASK=GLC  USED_N_YEARS=$(USED_N_YEARS) TAU_CCI=$(TAU_CCI) \
	APPLY_ABIOTIC_STATIC=$(APPLY_ABIOTIC_STATIC) \
	$(ENV_COMMON) $(R) $(WD)/07_apply_mask_0p05.R

# LUH overlays
.PHONY: apply-lai-cci-luhov
apply-lai-cci-luhov: cci-mask luh-overlap-mask abiotic-cci
	VAR=LAI MASK=CCI CCI_BAND=$(CCI_BAND) \
	APPLY_ABIOTIC_STATIC=$(APPLY_ABIOTIC_STATIC) APPLY_LUH_OVERLAP=TRUE \
	GRASS_SOURCE=$(GRASS_SOURCE) G_MIN=$(G_MIN) P_MIN=$(P_MIN) ALPHA=$(ALPHA) \
	LUH_AVG_START=$(LUH_AVG_START) LUH_AVG_END=$(LUH_AVG_END) TAU_CCI=$(TAU_CCI) \
	$(ENV_COMMON) $(R) $(WD)/07_apply_mask_0p05.R

.PHONY: apply-fpar-cci-luhov
apply-fpar-cci-luhov: cci-mask luh-overlap-mask abiotic-cci
	VAR=FPAR MASK=CCI CCI_BAND=$(CCI_BAND) \
	APPLY_ABIOTIC_STATIC=$(APPLY_ABIOTIC_STATIC) APPLY_LUH_OVERLAP=TRUE \
	GRASS_SOURCE=$(GRASS_SOURCE) G_MIN=$(G_MIN) P_MIN=$(P_MIN) ALPHA=$(ALPHA) \
	LUH_AVG_START=$(LUH_AVG_START) LUH_AVG_END=$(LUH_AVG_END) TAU_CCI=$(TAU_CCI) \
	$(ENV_COMMON) $(R) $(WD)/07_apply_mask_0p05.R

.PHONY: apply-lai-glc-luhov
apply-lai-glc-luhov: glc-mask luh-overlap-mask abiotic-cci
	VAR=LAI MASK=GLC USED_N_YEARS=$(USED_N_YEARS) \
	APPLY_ABIOTIC_STATIC=$(APPLY_ABIOTIC_STATIC) APPLY_LUH_OVERLAP=TRUE \
	GRASS_SOURCE=$(GRASS_SOURCE) G_MIN=$(G_MIN) P_MIN=$(P_MIN) ALPHA=$(ALPHA) \
	LUH_AVG_START=$(LUH_AVG_START) LUH_AVG_END=$(LUH_AVG_END) TAU_CCI=$(TAU_CCI) \
	$(ENV_COMMON) $(R) $(WD)/07_apply_mask_0p05.R

.PHONY: apply-fpar-glc-luhov
apply-fpar-glc-luhov: glc-mask luh-overlap-mask abiotic-cci
	VAR=FPAR MASK=GLC USED_N_YEARS=$(USED_N_YEARS) \
	APPLY_ABIOTIC_STATIC=$(APPLY_ABIOTIC_STATIC) APPLY_LUH_OVERLAP=TRUE \
	GRASS_SOURCE=$(GRASS_SOURCE) G_MIN=$(G_MIN) P_MIN=$(P_MIN) ALPHA=$(ALPHA) \
	LUH_AVG_START=$(LUH_AVG_START) LUH_AVG_END=$(LUH_AVG_END) TAU_CCI=$(TAU_CCI) \
	$(ENV_COMMON) $(R) $(WD)/07_apply_mask_0p05.R

# ==== aggregate masked 0.05° → 0.25° ====
.PHONY: agg-lai-cci
agg-lai-cci: apply-lai-cci
	VAR=LAI  MASK=CCI $(ENV_COMMON) $(R) $(WD)/08_agg_0p25.R

.PHONY: agg-fpar-cci
agg-fpar-cci: apply-fpar-cci
	VAR=FPAR MASK=CCI $(ENV_COMMON) $(R) $(WD)/08_agg_0p25.R

.PHONY: agg-lai-glc
agg-lai-glc: apply-lai-glc
	VAR=LAI  MASK=GLC $(ENV_COMMON) $(R) $(WD)/08_agg_0p25.R

.PHONY: agg-fpar-glc
agg-fpar-glc: apply-fpar-glc
	VAR=FPAR MASK=GLC $(ENV_COMMON) $(R) $(WD)/08_agg_0p25.R

.PHONY: agg-lai-cci-luhov
agg-lai-cci-luhov: apply-lai-cci-luhov
	VAR=LAI MASK=CCI $(ENV_COMMON) $(R) $(WD)/08_agg_0p25.R

.PHONY: agg-fpar-cci-luhov
agg-fpar-cci-luhov: apply-fpar-cci-luhov
	VAR=FPAR MASK=CCI $(ENV_COMMON) $(R) $(WD)/08_agg_0p25.R

.PHONY: agg-lai-glc-luhov
agg-lai-glc-luhov: apply-lai-glc-luhov
	VAR=LAI MASK=GLC $(ENV_COMMON) $(R) $(WD)/08_agg_0p25.R

.PHONY: agg-fpar-glc-luhov
agg-fpar-glc-luhov: apply-fpar-glc-luhov
	VAR=FPAR MASK=GLC $(ENV_COMMON) $(R) $(WD)/08_agg_0p25.R

# ==== high-level pipelines ====
.PHONY: lai-cci
lai-cci: config cci-mask abiotic-cci apply-lai-cci agg-lai-cci

.PHONY: fpar-cci
fpar-cci: config cci-mask abiotic-cci apply-fpar-cci agg-fpar-cci

.PHONY: lai-glc
lai-glc: config glc-mask abiotic-cci apply-lai-glc agg-lai-glc

.PHONY: fpar-glc
fpar-glc: config glc-mask abiotic-cci apply-fpar-glc agg-fpar-glc

.PHONY: lai-cci-luhov
lai-cci-luhov: config cci-mask luh-overlap-mask abiotic-cci apply-lai-cci-luhov agg-lai-cci-luhov

.PHONY: fpar-cci-luhov
fpar-cci-luhov: config cci-mask luh-overlap-mask abiotic-cci apply-fpar-cci-luhov agg-fpar-cci-luhov

.PHONY: lai-glc-luhov
lai-glc-luhov: config glc-mask luh-overlap-mask abiotic-cci apply-lai-glc-luhov agg-lai-glc-luhov

.PHONY: fpar-glc-luhov
fpar-glc-luhov: config glc-mask luh-overlap-mask abiotic-cci apply-fpar-glc-luhov agg-fpar-glc-luhov


# ==== quicklooks ====
.PHONY: ql
ql:
	REMAKE_QL=TRUE $(ENV_COMMON) $(R) $(WD)/02_cci_frac_0p05.R
	REMAKE_QL=TRUE $(ENV_COMMON) $(R) $(WD)/04_cci_mask_0p05.R
	REMAKE_QL=TRUE $(ENV_COMMON) $(R) $(WD)/03_glc_stack_0p05.R
	REMAKE_QL=TRUE $(ENV_COMMON) $(R) $(WD)/04_glc_mask_0p05.R

# ==== help ====
.PHONY: help
help:
	@echo "Available targets:"
	@echo "  all              : baseline pipeline (LAI/FPAR + CCI mask)"
	@echo "  all-cci-glc      : CCI + GLC pipelines (no LUH overlay)"
	@echo "  all-full         : CCI/GLC + LUH overlay"
	@echo "  precompute       : Run shared preprocessing (georef, CCI frac, GLC stack)"
	@echo "  lai-cci / fpar-cci / lai-glc / fpar-glc"
	@echo "  lai-cci-luhov / fpar-cci-luhov / lai-glc-luhov / fpar-glc-luhov"
	@echo "  luh-overlap-mask : create LUH overlap mask"
	@echo "  ql               : regenerate quicklooks"
	@echo "  RUN_TAG=$(RUN_TAG), COMMON_OUT=$(COMMON_OUT)"
	@echo "  CCI_BAND=$(CCI_BAND), TAU_CCI=$(TAU_CCI), K_CCI=$(K_CCI)"
	@echo "  GLC: USED_N_YEARS=$(USED_N_YEARS)"
	@echo "  Abiotic: $(TAU_WATER), $(TAU_ICE), $(TAU_BARE) from $(ABIOTIC_Y1) to $(ABIOTIC_Y2)"
	@echo "  LUH overlap: G_MIN=$(G_MIN), P_MIN=$(P_MIN), ALPHA=$(ALPHA), APPLY_LUH_OVERLAP=$(APPLY_LUH_OVERLAP)"
