# IgG N-Glycome Reconstruction

Reconstruction of IgG N-glycome profiles from the total blood plasma N-glycome (TPNG).

## Project Overview

This repository contains a pipeline for reconstructing immunoglobulin G (IgG) N-glycome profiles from the total plasma N-glycome (TPNG). The project includes data preparation, IgG peak prediction using linear and nonlinear models, reconstruction of IgG-specific GWAS summary statistics, meta-analysis, validation, and post-GWAS gene prioritization stages.

## Folders:

### 01_data_preparation:

Raw glycomics processing. Harmonization produced 36 plasma glycan peaks (PGP) and 21 IgG glycan peaks (IGP). Per cohort: probabilistic median quotient normalization, batch correction, log-transform, outlier removal (outside Q1 − 3×IQR / Q3 + 3×IQR), sex/age adjustment, and standardization (mean 0, SD 1).

### 02_peak_prediction:

IgG N-glycome reconstruction. Linear models: OLS baseline, elastic net (lasso/ridge), L2-regularized Huber, Theil–Sen. Nonlinear models: random forest, gradient boosted tree + linear regression, TabNet.

### 03_meta_analysis:

Meta-analysis scripts. 

### 04_validation:

Validation benchmark. 

### 05_postgwas:

Gene prioritization. 
