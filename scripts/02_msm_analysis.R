# ==============================================================================
# Script: 02_msm_analysis.R
# ==============================================================================

source(\"scripts/00_setup_and_bootstrap.R")

# 1. IPW Weight Creation Formula
# ------------------------------------------------------------------------------

compute_msm_weights_exact <- function(data) {
  
  # Ensure factors are appropriately classed
  data$CF <- as.numeric(as.character(data$CF))
  data$MV <- as.numeric(as.character(data$MV))
  
  # A. Base Models (Visit == 0)
  df_baseline <- filter(data, visit == 0)
  mod_num0 <- glm(polypharmacy_med ~ 1, data = df_baseline, family = binomial())
  
  den0_formula <- formula(
    polypharmacy_med ~ gender + educ3cat + agerand + CF + MV +
      BMI + smoke2 + hxhtn + diabetes + depression + cvd + cancer + 
      currhealth + arrhythmia + thrombosis + Cirrhosis_liver + parkinson +
      diabetes*gender + diabetes*hxhtn + depression*gender +
      depression*currhealth + depression*diabetes + gender*agerand +
      gender*currhealth + depression*smoke2 + hxhtn*depression +
      hxhtn*gender + hxhtn*agerand + hxhtn*BMI + hxhtn*smoke2 +
      hxhtn*educ3cat + hxhtn*arrhythmia
  )
  mod_den0 <- glm(den0_formula, data = df_baseline, family = binomial())
  
  data <- data %>%
    mutate(
      pnum0   = if_else(visit == 0, predict(mod_num0, newdata = ., type = "response"), NA_real_),
      pdenom0 = if_else(visit == 0, predict(mod_den0, newdata = ., type = "response"), NA_real_)
    )
  
  # B. Time-Varying Models (Visit >= 1)
  nFit <- glm(polypharmacy_med ~ visit + visit2 + poly_lag1, data = data, family = binomial())
  
  # Exact dFit matching baseline (_b), time-varying, and lagged interactions
  dFit_formula <- formula(
    polypharmacy_med ~ visit + visit2 + 
      gender + educ3cat + agerand + CF + MV + 
      BMI_b + smoke2_b + hxhtn_b + diabetes_b + depression_b + cvd_b + cancer_b + 
      currhealth_b + arrhythmia_b + thrombosis_b + kidney_failure_b + Cirrhosis_liver_b + parkinson_b +
      BMI + smoke2 + hxhtn + diabetes + depression + cvd + cancer + currhealth + arrhythmia +
      thrombosis + kidney_failure + Cirrhosis_liver + parkinson +
      BMI_lag + smoke2_lag + hxhtn_lag + diabetes_lag + depression_lag +
      cvd_lag + cancer_lag + currhealth_lag + arrhythmia_lag +
      thrombosis_lag + kidney_failure_lag + Cirrhosis_liver_lag + parkinson_lag +
      gender*BMI_lag + diabetes_lag*agerand + diabetes_lag*BMI_lag + diabetes_lag*educ3cat +
      depression_lag*hxhtn_lag + currhealth_lag*depression_lag + poly_lag1
  )
  
  dFit <- glm(dFit_formula, data = data, family = binomial())
  
  data <- data %>%
    mutate(
      pnum   = if_else(visit != 0, predict(nFit, newdata = ., type = "response"), NA_real_),
      pdenom = if_else(visit != 0, predict(dFit, newdata = ., type = "response"), NA_real_)
    )
  
  # C. Cumulative Weights and Truncation
  data <- data %>%
    mutate(
      numCont = if_else(visit == 0, if_else(polypharmacy_med == 1, pnum0, 1 - pnum0), if_else(polypharmacy_med == 1, pnum, 1 - pnum)),
      denCont = if_else(visit == 0, if_else(polypharmacy_med == 1, pdenom0, 1 - pdenom0), if_else(polypharmacy_med == 1, pdenom, 1 - pdenom))
    ) %>%
    group_by(id) %>%
    mutate(cum_num = cumprod(numCont), cum_den = cumprod(denCont), stabw = cum_num / cum_den) %>%
    ungroup()
  
  t0 <- quantile(data$stabw[data$visit == 0], 1, na.rm = TRUE)
  t1 <- quantile(data$stabw[data$visit >= 1], 1, na.rm = TRUE)
  
  data %>% mutate(stabw_t = case_when(visit == 0 & stabw > t0 ~ t0, visit >= 1 & stabw > t1 ~ t1, TRUE ~ stabw))
}

# 2. Pooled Logistic Model & Bootstrap CI Calculation
# ------------------------------------------------------------------------------
# The Marginal Structural Model (MSM) utilizes the weights (stabw_t) computed above.
# The core outcome formula evaluated leverages visit-time interactions:
# Outcome ~ polypharmacy_med + visit + I(visit^2) + (polypharmacy_med * visit) + (polypharmacy_med * visit^2)
#
# Note: The extensive bootstrapping loop (B=500 iterations) required to calculate the
# precise 95% Confidence Intervals for Cumulative Incidence Ratios (CIR) is defined
# as a core function `run_msm_bootstrap_CIR` in scripts/00_setup_and_bootstrap.R.

# --- Execution ---


# 1. Define the outcome formula (e.g., All-Cause Death)
formula_allcause <- "DEATH ~ polypharmacy_med * (visit + I(visit^2))"

# 2. Execute the Bootstrap
msm_results_allcause <- run_msm_bootstrap_CIR(
  data = df_analysis, # Replace 'df_analysis' with your dataset
  formula_str = formula_allcause,
  weight_func = compute_msm_weights_exact,
  weight_col = "stabw_t",
  B = 500
)
print(msm_results_allcause)

cat("IPW creation logic defined. MSM Bootstrapping functions are available via 00_setup_and_bootstrap.R\n")
