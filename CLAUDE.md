# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based research project analyzing arbuscular mycorrhizal fungi (AMF) community assembly in coffee-growing regions of Uganda (ROBUST project). It integrates 18S rRNA metabarcoding data, remote sensing (NDVI from Sentinel-3 via Planetary Computer), and spatial covariates to build Joint Species Distribution Models (JSDMs).

## Running the Analysis

There is no build system — scripts are run sequentially in RStudio by knitting `.Rmd` files or sourcing `.R` files:

```r
# Run a specific step
rmarkdown::render("code/3_rarefaction_AMF.Rmd")

# Or source R scripts directly
source("code/5_jsdm_genus_models.R")
```

The pipeline must be run in order (1 → 9) as each step produces outputs consumed by the next.

## Pipeline Architecture

| Script | Role |
|--------|------|
| `1_extract_raster_data.Rmd` | Fetches NDVI rasters from Planetary Computer STAC API (Sentinel-3 Synergy, Uganda region) |
| `2_raster_formatting.Rmd` | Preprocesses and formats raster data |
| `3_rarefaction_AMF.Rmd` | Filters Glomeromycota/Mucoromycotina OTUs from 18S MiSeq data; rarefaction |
| `4_setup_data_for_model.Rmd` | Merges spatial covariates with sample metadata; UTM → WGS84 conversion |
| `5_jsdm_genus_models.R` | Fits all 32 additive Hmsc JSDMs (Gaussian + probit) with phylo tree and spatial + iid random levels |
| `6_convergence_diagnostics.Rmd` | Convergence checks (PSRF, WAIC) and DHARMa residual diagnostics |
| `7_cross_validation_best_models.Rmd` | Cross-validation of best WAIC models (RMSE for Gaussian, AUC for probit) |
| `8_jsdm_results.Rmd` | Combined Gaussian + Probit model results, merged effect/phylo-signal figures |
| `9_descriptives_figure.Rmd` | Descriptive statistics and maps |

### Deprecated

| Script | Reason |
|--------|--------|
| `deprecated_6_jsdm_genus_final.Rmd` | Replaced by script 8; loads old model filenames |
| `deprecated_5_timing_test.R` | Benchmark utility, not part of the pipeline |
| `deprecated_8_jsdm_gauss_results.Rmd` | Merged into `8_jsdm_results.Rmd` |
| `deprecated_9_jsdm_bin_results.Rmd` | Merged into `8_jsdm_results.Rmd` |

## Key Data Paths

- **Raw OTU table**: `data/18S_AMf_Eukaryome/Robust_MiSeq_18S_AMf_20241120_125_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2`
- **Raster data**: `data/raster/`
- **Model outputs**: `outputs/models/`
- **Processed abundance tables**: `outputs/abundance_tables/`

Note: Most data files are git-ignored (`.csv`, `.txt`, `.tif`, `.png`, etc.).

## Core R Packages

- **Hmsc** — Hierarchical Modelling of Species Communities (JSDM fitting)
- **rstac** — Planetary Computer STAC API access for raster fetching
- **terra / stars / raster** — Raster spatial operations
- **sf / sp** — Vector spatial data and coordinate transforms
- **vegan** — Rarefaction and community ecology
- **ape / ggtree** — Phylogenetic tree handling and visualization
- **tidyverse** — Data wrangling and plotting

## Custom Code

- `code/ggheatmap.custom.R` — Modified `ggheatmap` function that overlays heatmaps on phylogenetic trees (extends `ggtree`). Source this before scripts that produce phylo-heatmap figures.

## Code Style

- RStudio project uses **2-space tab indentation**, UTF-8 encoding.
- Scripts use `tidyverse` conventions (pipes, `dplyr` verbs).
