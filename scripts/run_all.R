# ==============================================================================
# Script: run_all.R
# Author: COSMOS Polypharmacy Research Team
# ==============================================================================

# Time the full run
start_time <- Sys.time()
cat(sprintf("Starting execution at %s\n", start_time))
cat("==========================================================================\n")

# Execute pipelines sequentially
cat("1/8: Loading setup functions...\n")
source(\"scripts/00_setup_and_bootstrap.R")

cat("2/8: Prepping Baseline Cohort & Processing Outcomes...\n")
source(\"scripts/01_data_cleaning_and_outcomes.R")

cat("3/8: Unfurling Follow-Up Iterations (Longitudinal Processing)...\n")
source(\"scripts/02_longitudinal_processing.R")

cat("4/8: Analyzing Missingness & Imputation...\n")
source(\"scripts/03_imputation.R")

cat("5/8: Running Core Cox Proportional Hazard Models...\n")
source(\"scripts/01_cox_regression.R")

cat("6/8: Evaluating Dynamic IPW with Marginal Structural Models (MSM)...\n")
source(\"scripts/02_msm_analysis.R")

cat("7/8: Expanding to Subgroups and Sensitivity Modules...\n")
source(\"scripts/03_subgroup_analyses.R")
source(\"scripts/04_sensitivity_analyses.R")

cat("8/8: Assembling Tables, Outputs, and High-Res Plots...\n")
source(\"scripts/05_tables_and_figures.R")

# Done
end_time <- Sys.time()
diff_time <- round(difftime(end_time, start_time, units = "mins"), 2)
cat("==========================================================================\n")
cat(sprintf("SUCCESS: End-to-end analytical pipeline completed flawlessly in %s minutes.\n", diff_time))
