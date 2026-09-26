# Bachelor thesis

**Robustness of DMI-Informed Genome-Scale Metabolic Modeling to Transcriptomic Reconstruction Choices**

This repository contains the code, data, models, results, figures, and tables used in my bachelor thesis.

## Contents

- `data/` - transcriptomic and DMI input data
- `scripts/` - analysis code and notebooks
- `models/` - metabolic models
- `results/` - analysis results
- `figures/` - figures used in the thesis
- `tables/` - tables and source data
- `config/` - model settings
- `environment/` - software environment information
- `tests/` - basic checks

## Data

The transcriptomic datasets used in the thesis are GSE17576 and GSE305719.

The original GEO files are in `data/raw_public/`.

Processed expression data and sample information are in `data/processed_transcriptomics/`.

DMI input data are in `data/dmi/`.

## Analysis

The repository includes code for:

- transcriptomic preprocessing
- gene-confidence assignment
- CORDA model reconstruction at P5, P10, P15, P20, P25, P30, P35, P40 and P45
- DMI data processing and model fitting
- integration of DMI-derived MR_glc and V_ox values into the metabolic models
- model setup and feasibility testing
- flux balance analysis
- flux variability analysis
- flux sampling
- reaction and pathway comparisons
- comparison of HFD and control conditions
- comparison across reconstruction thresholds
- comparison between GSE17576 and GSE305719
- statistical analysis
- generation of figures and tables

The DMI analysis code is in `scripts/dmi/`.

Transcriptomic processing code is in `scripts/transcriptomics/`.

## Public transcriptomic data

GSE17576 and GSE305719 are available from the Gene Expression Omnibus.

The DMI-derived data originate from experimental measurements described by Pang et al. (2026).
