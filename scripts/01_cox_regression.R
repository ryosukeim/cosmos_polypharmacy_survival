# ==============================================================================
# Script: 01_cox_regression.R
# ==============================================================================

source(\"scripts/00_setup_and_bootstrap.R")

df_analysis <- readRDS("data/processed/analysis_dataset.rds")

# 1. Prepare Baseline specific for Cox
# ------------------------------------------------------------------------------
d_base <- df_analysis %>% filter(month == 0)

# 2. Cox Regression Variables
# ------------------------------------------------------------------------------
# Time scales: t2death/12, etc based on original outcome loops
run_cox_pipeline <- function(d, survival_obj, exposure_var) {
  covars <- c("agerand", "gender", "cvd", "educ3cat", "BMI", "smoke", 
              "hxhtn", "diabetes", "arrhythmia", "thrombosis", "Cirrhosis_liver", 
              "parkinson", "kidney_failure", "depression", "cancer", "MV", "CF")
  
  form_crude <- as.formula(paste(survival_obj, "~", exposure_var))
  form_adj   <- as.formula(paste(survival_obj, "~", exposure_var, "+", paste(covars, collapse = "+")))
  
  crude_fit <- coxph(form_crude, data = d)
  adj_fit   <- coxph(form_adj, data = d)
  
  list(crude = crude_fit, adjusted = adj_fit)
}

# 3. Execution
# ------------------------------------------------------------------------------
cat("Running Cox Regressions...\n")

# All-Cause Death
cox_allcause <- run_cox_pipeline(d_base, "Surv(t2death/12, DEATH)", "polypharmacy_med")
cox_allcause_qt <- run_cox_pipeline(d_base, "Surv(t2death/12, DEATH)", "factor(num_med_quartile)") # Assuming computed

# CVD Death
cox_cvd <- run_cox_pipeline(d_base, "Surv(t2cvddth/12, cvddth)", "polypharmacy_med")
# Cancer Death
cox_cancer <- run_cox_pipeline(d_base, "Surv(t2cancdth/12, cancdth)", "polypharmacy_med")
# Other Death (Non-CVD / Non-Cancer)
cox_other <- run_cox_pipeline(d_base, "Surv(t2mistrrevcvddthoth/12, mistrrevcvddthoth)", "polypharmacy_med")

# Results will be collected in 08_tables
cat("Cox Proportional Hazards array configured.\n")
