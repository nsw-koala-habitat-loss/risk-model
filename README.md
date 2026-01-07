# Modelling deforestation risk and its implications on koala habitats

Author: Yi Fei Chung & Jonathan R. Rhodes
Last updated: 05 January 2026

Code to reproduce the analyses and figures in the manuscript (DOI:).

## Overview
1. All covariates used in this analysis were processed in the "nsw-koala-habitat-loss/risk-model-covariates" repository.
2. These covariates were then preprocessed to aggregated to cadastral land parcel units in [preprocess script](pre_processing.R)
3. The main analysis is performed in the [Model fitting](model_fitting.R) script.
4. Figures and tables are generated in the [plotting ](plotting.R) script.

## INLA posterior sampling
Predictions step involves drawing 50 000 samples from the posterior distribution approximated by the final selected model. Memory requirements for this step can be >1TB for certain models (e.g. infrastructure clearing model for Central Coast region). Bunya HPC was used to perform this step.

As the computing resources required are different for each model, a list of rscripts and scheduling scripts (SLURM) were used for executing the prediction step. See [Generate R script](Generate_R_script.R) and  [Generate SLURM script](Generate_SLURM_Script.R) for estimating HPC resource and generating SLURM script for HPC job scheduling. To avoid redundancy, actual scripts were not included.
