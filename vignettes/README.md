# **Vignettes**

This directory contains extended documentation for the natural-vegetation LAI/FPAR processing workflow.
The vignettes provide a narrative overview of the masking, aggregation, and evaluation pipeline, complementing the step-wise R scripts under `R/`.

## **Contents**

* **`vignette.Rmd`**
  A full workflow description covering:

  * configuration and environment setup
  * mask construction (CCI, GLC, LUH, abiotic)
  * application of masks to monthly LAI/FPAR
  * area-weighted aggregation to 0.25°
  * sensitivity analysis and evaluation
  * visual exploration of outputs

* **`README.md`**
  This file. Summarizes the purpose and structure of the vignette directory.

## **Usage**

The R Markdown vignette is designed to be run interactively or rendered via:

```r
rmarkdown::render("vignette.Rmd")
```

Scripts referenced in the vignette are fully non-interactive (HPC-ready) and can be invoked directly using `Rscript` or via the Makefile.

## **Scope**

The vignettes are not execution units.
They serve as:

* an overview of the conceptual workflow
* documentation for reproducibility
* guidance for new users or collaborators
* a reference for mask semantics, aggregation rules, and evaluation practices

All computations remain in the main pipeline (under `R/`), ensuring separation between documentation and production code.

---

