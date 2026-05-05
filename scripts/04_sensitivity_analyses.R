# ==============================================================================
# Script: 04_sensitivity_analyses.R
# ==============================================================================

source(\"scripts/00_setup_and_bootstrap.R")

cat("Configuring Sensitivity Paradigms...\n")

# 1. Comorbidity == 1 Restriction
# ------------------------------------------------------------------------------
get_comorbidity_cohort <- function(df_analysis) {
  com_yes_ids <- df_analysis %>% 
    filter(month == 0) %>% 
    filter(hxhtn == 1 | diabetes == 1 | depression == 1 | cvd == 1 | cancer == 1 | 
             arrhythmia == 1 | thrombosis == 1 | kidney_failure == 1 | 
             Cirrhosis_liver == 1 | parkinson == 1) %>% 
    pull(id)
  
  df_com_yes <- df_analysis %>% filter(id %in% com_yes_ids)
  return(df_com_yes)
}

# 2. Excluding Comorbidities from IPW Model
# ------------------------------------------------------------------------------
# Drops cardiovascular and neurological parameters from denominator formula
get_den0_excluding_comorb <- function(df_baseline) {
  glm(polypharmacy_med ~ gender + educ3cat + agerand + CF + MV + BMI + smoke2 + currhealth,
      data = df_baseline, family = binomial())
}

# 3. Incorporating Supplements (Exposure: polypharmacy_with_supplements instead of med)
# ------------------------------------------------------------------------------
# Shifts primary weight creation models to polypharmacy_with_supplements and polypharmacy_with_supplements_lag
get_den0_supplements <- function(df_baseline) {
  glm(polypharmacy_with_supplements ~ gender + educ3cat + agerand + CF + MV +
        BMI + smoke2 + hxhtn + diabetes + depression + cvd + cancer + 
        currhealth + arrhythmia + thrombosis + Cirrhosis_liver + parkinson +
        diabetes*gender + diabetes*hxhtn + gender*diabetes + depression*gender +
        depression*currhealth + depression*diabetes + gender*agerand +
        gender*currhealth + depression*smoke2,
      data = df_baseline, family = binomial())
}

# 4. Potentially Inappropriate Medications (PIM) MULTINOMIAL Weights
# ------------------------------------------------------------------------------
get_pim_multinomial_weights <- function(df_analysis) {
  
  # Ensure polypharmacy_type is factored
  df_analysis <- df_analysis %>%
    mutate(
      polypharmacy_type = factor(polypharmacy_type, levels = c("0_NonPoly", "1_PIM_free", "2_PIM_containing")),
      poly_type_chr = as.character(polypharmacy_type)
    )
  
  d0 <- filter(df_analysis, visit == 0)
  
  numFit0 <- nnet::multinom(polypharmacy_type ~ 1, data = d0, trace = FALSE)
  
  denFit0 <- nnet::multinom(polypharmacy_type ~ gender + educ3cat + agerand + CF + MV +
                              BMI + smoke2 + hxhtn + diabetes + depression + cvd + cancer + 
                              currhealth + arrhythmia + thrombosis + Cirrhosis_liver + parkinson +
                              diabetes*gender + diabetes*hxhtn + gender*diabetes + depression*gender +
                              depression*currhealth + depression*diabetes + gender*agerand +
                              gender*currhealth + cvd*gender + depression*smoke2,
                            data = d0, trace = FALSE)
  
  pnum0 <- predict(numFit0, newdata = d0, type = "probs")
  pden0 <- predict(denFit0, newdata = d0, type = "probs")
  
  # Extracts multinomial propensity scoring for exact observed group
  extract_prob <- function(prob_vec, observed_level) prob_vec[observed_level]
  
  d0$numCont_pim <- mapply(function(i, obs) extract_prob(pnum0[i, ], obs), seq_len(nrow(d0)), d0$poly_type_chr)
  d0$denCont_pim <- mapply(function(i, obs) extract_prob(pden0[i, ], obs), seq_len(nrow(d0)), d0$poly_type_chr)
  
  d0 <- d0 %>%
    group_by(id) %>%
    mutate(
      k1_0_pim = cumprod(numCont_pim),
      k1_w_pim = cumprod(denCont_pim),
      stabw_pim = k1_0_pim / k1_w_pim
    ) %>%
    ungroup()
  
  threshold0 <- quantile(d0$stabw_pim, 1, na.rm = TRUE)
  d0 <- d0 %>% mutate(stabw_pim_t = ifelse(stabw_pim > threshold0, threshold0, stabw_pim))
  
  # Note: A similar multinomial regression must cover d_rest (visit > 0)
  return(d0)
}



# 5. Exclude Aspirin & NSAIDs
# ------------------------------------------------------------------------------
# Converts primary exposure object iteratively to polypharmacy_med_wo_otc and its lag
get_den0_exclude_aspirin <- function(df_baseline) {
  glm(polypharmacy_med_wo_otc ~ gender + educ3cat + agerand + CF + MV +
        BMI + smoke2 + hxhtn + diabetes + depression + cvd + cancer + 
        currhealth + arrhythmia + thrombosis + Cirrhosis_liver + parkinson +
        diabetes*gender + diabetes*hxhtn + gender*diabetes + depression*gender +
        depression*currhealth + depression*diabetes + gender*agerand +
        gender*currhealth + depression*smoke2 +
        hxhtn*depression + hxhtn*gender + hxhtn*agerand + hxhtn*BMI + hxhtn*smoke2 +
        hxhtn*educ3cat + hxhtn*arrhythmia,
      data = df_baseline, family = binomial())
}

# 6. Additional Covariates (Sleep, Income, PU)
# ------------------------------------------------------------------------------
# Integrates `income`, plus `PU`, `sleep_less` and their corresponding lags into the time-varying matrix denominator
get_den0_additional_covariates <- function(df_baseline) {
  glm(polypharmacy_med ~ gender + educ3cat + agerand + CF + MV + income +
        BMI + smoke2 + hxhtn + diabetes + depression + cvd + cancer + 
        currhealth + arrhythmia + thrombosis + Cirrhosis_liver + parkinson +
        PU + sleep_less +
        diabetes*hxhtn + gender*diabetes + depression*gender +
        depression*currhealth + depression*diabetes + gender*agerand +
        gender*currhealth + depression*smoke2,
      data = df_baseline, family = binomial())
}

# 7. Truncating weights at the 99.5th percentile
# ------------------------------------------------------------------------------
# Limits the standard calculated inverse probability weight vector at 0.995
apply_truncation_995 <- function(df_analysis) {
  threshold0 <- quantile(df_analysis$stabw[df_analysis$visit == 0], 0.995, na.rm = TRUE)
  threshold1 <- quantile(df_analysis$stabw[df_analysis$visit >= 1], 0.995, na.rm = TRUE)
  
  df_analysis %>%
    mutate(stabw_t995 = case_when(
      visit == 0  & stabw > threshold0 ~ threshold0,
      visit >= 1  & stabw > threshold1 ~ threshold1,
      TRUE ~ stabw
    ))
}

# 8. Adjusting for loss to follow-up using Censoring IPCW
# ------------------------------------------------------------------------------
get_censoring_ipcw_weights <- function(df_analysis, data_death) {
  
  # 1. Join death times and define the 'uncensored' indicator
  # Use orig_id if it exists (e.g., from bootstrap resampling), otherwise use id
  join_col <- if ("orig_id" %in% names(df_analysis)) "orig_id" else "id"
  
  d_uncen <- df_analysis %>%
    left_join(data_death %>% select(id, t2death) %>% rename(!!join_col := id), by = join_col) %>%
    mutate(t2death_month = t2death / 30.44) %>%
    arrange(id, month) %>%
    group_by(id) %>%
    mutate(
      # 'uncensored0' remains 1 while still under follow-up
      uncensored0 = if_else(!is.na(month) & month <= t2death_month, 1L, 0L),
      uncensored  = cumprod(uncensored0)
    ) %>%
    ungroup()
    
  # 2. Fit censoring models for follow-up visits (visit >= 1)
  df_fit <- d_uncen %>% filter(visit >= 1)
  
  num_fit_uncensored <- glm(
    uncensored ~ visit + I(visit^2) + poly_lag1,
    data   = df_fit,
    family = binomial()
  )
  
  den_fit_uncensored <- glm(
    uncensored ~ visit + I(visit^2) + gender + educ3cat + agerand + CF + MV +
      BMI_b + smoke2_b + hxhtn_b + diabetes_b + depression_b + cvd_b + cancer_b + 
      currhealth_b + arrhythmia_b + thrombosis_b + kidney_failure_b + 
      Cirrhosis_liver_b + parkinson_b + BMI + smoke2 + hxhtn + diabetes + 
      depression + cvd + cancer + currhealth + arrhythmia + thrombosis + 
      kidney_failure + Cirrhosis_liver + parkinson + BMI_lag + smoke2_lag + 
      hxhtn_lag + diabetes_lag + depression_lag + cvd_lag + cancer_lag + 
      currhealth_lag + arrhythmia_lag + thrombosis_lag + kidney_failure_lag + 
      Cirrhosis_liver_lag + parkinson_lag + poly_lag1,
    data   = df_fit,
    family = binomial()
  )
  
  # 3. Predict per-visit probabilities
  d_uncen <- d_uncen %>%
    mutate(
      pnum0_uncensored = if_else(visit == 0, 1, NA_real_),
      pden0_uncensored = if_else(visit == 0, 1, NA_real_),
      pnum_uncensored  = if_else(visit != 0, predict(num_fit_uncensored, newdata = ., type = "response"), NA_real_),
      pden_uncensored  = if_else(visit != 0, predict(den_fit_uncensored, newdata = ., type = "response"), NA_real_)
    )
    
  # 4. Compute per-visit weight contributions
  d_uncen <- d_uncen %>%
    mutate(
      numCont_uncensored = case_when(
        visit == 0 ~ if_else(uncensored == 1, pnum0_uncensored, 1 - pnum0_uncensored),
        TRUE       ~ if_else(uncensored == 1, pnum_uncensored,  1 - pnum_uncensored)
      ),
      denCont_uncensored = case_when(
        visit == 0 ~ if_else(uncensored == 1, pden0_uncensored, 1 - pden0_uncensored),
        TRUE       ~ if_else(uncensored == 1, pden_uncensored,  1 - pden_uncensored)
      )
    )
    
  # 5. Accumulate and derive IPCW, truncate at 99.5th percentile
  d_uncen <- d_uncen %>%
    group_by(id) %>%
    mutate(
      cum_num_uncensored = cumprod(numCont_uncensored),
      cum_den_uncensored = cumprod(denCont_uncensored),
      ipcw_uncensored    = cum_num_uncensored / cum_den_uncensored
    ) %>%
    ungroup()
    
  thr_uncensored <- quantile(d_uncen$ipcw_uncensored, probs = 0.995, na.rm = TRUE)
  
  d_uncen <- d_uncen %>%
    mutate(
      ipcw_uncensored_t = pmin(ipcw_uncensored, thr_uncensored),
      pp_wt_uncensored  = stabw_t * ipcw_uncensored_t
    )
    
  return(d_uncen)
}

cat("Sensitivity definitions configured.\n")
