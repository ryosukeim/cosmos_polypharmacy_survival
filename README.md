Polypharmacy and Mortality in Older Adults: Time-Varying Exposure Analysis of the COSMOS Study

Overview
This repository contains the analytical code for the paper: "Polypharmacy and mortality in older adults: time-varying exposure analysis of the COSMOS study".

The study investigates the longitudinal association between persistent polypharmacy (concurrent use of >= 5 medications) and all-cause mortality among community-dwelling older adults, utilizing inverse-probability weighting (IPW) to adjust for time-varying confounding.


Repository Structure
- scripts/00_setup_and_bootstrap.R: Custom R functions and package dependencies.
- scripts/01_cox_regression.R: Standard unadjusted and adjusted Cox proportional hazard models.
- scripts/02_msm_analysis.R: Marginal Structural Models for time-varying confounding.
- scripts/03_subgroup_analyses.R: Effect modification across subgroups.
- scripts/04_sensitivity_analyses.R: Sensitivity analysis.
- scripts/05_tables_and_figures.R: Tables and figures.
- scripts/run_all.R: Master execution script.

Requirements
R version 4.4.2
Packages: tidyverse, haven, survival, survey, mice, nnet, tableone
