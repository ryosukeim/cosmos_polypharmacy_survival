# ==============================================================================
# Script: 05_tables_and_figures.R
# ==============================================================================

source(\"scripts/00_setup_and_bootstrap.R")

cat("Executing Result Formatting and Table Generation...\n")

# 1. Unweighted and IPW-Weighted Table 1
# ------------------------------------------------------------------------------
generate_table_one_weighted <- function(df_analysis, exposure_var = "polypharmacy_med") {
  
  vars <- c("gender", "agerand", "BMI", "educ3cat", "smoke2", "currhealth",  
            "hxhtn", "cvd", "diabetes",  "arrhythmia", "thrombosis", "kidney_failure", 
            "Cirrhosis_liver", "parkinson","depression", "cancer", 
            "CF", "MV")
  
  factorVars <- c("gender", "educ3cat", "smoke2", 
                  "hxhtn", "diabetes", "depression", "arrhythmia",
                  "cvd", "cancer","CF", "MV",
                  "thrombosis", "kidney_failure", "Cirrhosis_liver", "parkinson")
  
  # Ensure formats match
  df_analysis[factorVars] <- lapply(df_analysis[factorVars], factor)
  
  # Survey Design object for IPW adjustment (using stabw_t)
  if("stabw_t" %in% names(df_analysis)) {
    design0 <- subset(svydesign(ids = ~1, data = df_analysis, weights = ~stabw_t), month == 0)
    
    tab1_weighted0 <- svyCreateTableOne(vars = vars, strata = exposure_var, 
                                        data = design0, factorVars = factorVars)
    
    print(tab1_weighted0, smd = TRUE, explain = FALSE, varLabels = TRUE) %>% as.data.frame()
  } else {
    stop("stabilized IPW weights `stabw_t` not found in dataset for weighting.")
  }
}

# 2. Extract Baseline Non-Polypharmacy Cumulative Incidence
# ------------------------------------------------------------------------------
# Calculates `np_cuminc` (base absolute risk) needed to anchor relative CIRs
extract_np_cuminc <- function(df_analysis, formula_str, weights_col = "stabw_t") {
  plr_all <- glm(as.formula(formula_str),
                 data    = df_analysis,
                 weights = df_analysis[[weights_col]],
                 family  = quasibinomial())
  
  base <- df_analysis %>% filter(visit == 0)
  n    <- nrow(base)
  reps <- 0:4
  
  df0 <- base[rep(seq_len(n), each = length(reps)), ] %>%
    mutate(
      visit            = rep(reps, times = n),
      visit2           = visit^2,
      polypharmacy_med = 0,
      poly_med0visit   = 0,
      poly_med0visit2  = 0
    )
    
  one_minus_p <- 1 - predict(plr_all, newdata = df0, type = "response")
  df0$s <- ave(one_minus_p, df0$id, FUN = cumprod)
  
  nonpoly_inc <- 1 - tapply(df0$s, df0$visit, mean)
  return(nonpoly_inc[-1]) # return visits 1-4
}

# 3. Generating Weighted Cumulative Incidence Curves 
# ------------------------------------------------------------------------------
plot_cumulative_incidence_curve <- function(CI_tbl, np_cuminc, outcome_label = "all-cause mortality") {
  cir_pt <- CI_tbl$CIR

  poly_inc_pt <- (np_cuminc * cir_pt)
  np_inc_pt   <- np_cuminc
  
  plot_df <- tibble(
    visit = rep(0:4, 2),
    arm   = rep(c("Scenario B: Sustained No-polypharmacy", "Scenario A: Sustained Polypharmacy"), each = 5),
    est   = c(0, np_inc_pt, 0, poly_inc_pt)
  )
  
  ggplot(plot_df, aes(x = visit, y = est, colour = arm, fill = arm)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    scale_x_continuous(breaks = 0:4) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_colour_manual(values = c("Scenario B: Sustained No-polypharmacy" = "#0072B5FF", "Scenario A: Sustained Polypharmacy" = "#BC3C29FF")) +
    scale_fill_manual(values = c("Scenario B: Sustained No-polypharmacy" = "#0072B5FF", "Scenario A: Sustained Polypharmacy" = "#BC3C29FF")) +
    labs(x = "Time (years)", y = paste("Cumulative Incidence of", outcome_label), colour = NULL, fill = NULL) +
    theme_bw(base_size = 14) + theme(legend.position = "bottom")
}

# 4. Generating Quartile KM Curves (Cumulative Incidence format)
# ------------------------------------------------------------------------------
# Uses survminer to plot baseline quartile comparisons accurately natively
plot_quartile_survival_curve <- function(surv_fit, data, max_y = 0.1) {
  # Requires library(survminer) and library(survival)
  p <- survminer::ggsurvplot(
    surv_fit,
    data = data,
    fun = "event",
    risk.table = TRUE,
    pval = FALSE,
    legend.title = "",
    legend.labs = c("Q1 (0-1)", "Q2 (2)", "Q3 (3-4)", "Q4 (>=5)"),
    xlab = "Time (years)", 
    ylab = "Cumulative incidence of all-cause death (%)",
    palette = c("#0072B5FF", "#20854EFF", "#E18727FF", "#BC3C29FF"),
    censor = FALSE,
    xlim = c(0, 4),
    break.time.by = 1
  )
  
  p$plot <- p$plot + scale_y_continuous(
    labels = scales::label_percent(accuracy = 0.1),
    limits = c(0, max_y),
    breaks = seq(0, max_y, by = 0.02),
    expand = c(0, 0)
  ) + theme(legend.text = element_text(size = 14))
  
  return(p)
}

# 5. Generating Baseline Binary KM Curves (Cumulative Incidence format)
# ------------------------------------------------------------------------------
plot_baseline_binary_survival_curve <- function(surv_fit, data, max_y = 0.1) {
  # Requires library(survminer) and library(survival)
  p <- survminer::ggsurvplot(
    surv_fit,
    data = data,
    fun = "event",
    risk.table = TRUE,
    pval = FALSE,
    legend.title = "",
    legend.labs = c("Scenario B: Sustained No-polypharmacy", "Scenario A: Sustained Polypharmacy"),
    xlab = "Time (years)", 
    ylab = "Cumulative incidence of all-cause death (%)",
    palette = c("#0072B5FF", "#BC3C29FF"),
    censor = FALSE,
    xlim = c(0, 4),
    break.time.by = 1
  )
  
  p$plot <- p$plot + scale_y_continuous(
    labels = scales::label_percent(accuracy = 0.1),
    limits = c(0, max_y),
    breaks = seq(0, max_y, by = 0.02),
    expand = c(0, 0)
  ) + theme(legend.text = element_text(size = 14))
  
  return(p)
}

cat("Table 1 and Graphical weighting paradigms configured using package `survey` and `tableone`.\n")
