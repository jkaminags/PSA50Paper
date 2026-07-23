library(tidyverse)
# int.simpson2() (AUC estimation) is called via fda.usc::. Attaching the whole
# package masks dplyr::select, which breaks the test functions below.
source("polynomial_PSA.R")

# Generates a simulated PSA database
generate_PSA_sim <- function(
    baseline_mean_A, baseline_sd_A,
    baseline_mean_B, baseline_sd_B,
    response_prop_A, response_prop_B,
    correlation_sd
) {
  groups <- rep(c(rep("A", 30), rep("B", 30)), 8)
  patient_id <- rep(seq(1:60), 8)
  time <- c(rep(0, 60), rep(2, 60), rep(3, 60), rep(4, 60),
            rep(6, 60), rep(8, 60), rep(9, 60), rep(12, 60))
  
  # baselines floored at 0.1 ng/mL for the same reason (see generate_PSA_FUP)
  baseline_A <- pmax(rnorm(30, baseline_mean_A, baseline_sd_A), 0.1)
  baseline_B <- pmax(rnorm(30, baseline_mean_B, baseline_sd_B), 0.1)
  
  # generate grp A follow-up PSA at 3, 6, 9, and 12 months
  # add in measurement timepoints to reflect CASCARA data
  
  fup_A_2 <- generate_PSA_FUP(baseline_mean_A, baseline_A,
                              response_prop_A, correlation_sd, 2)
    
  fup_A_3 <- generate_PSA_FUP(baseline_mean_A, baseline_A,
                            response_prop_A, correlation_sd, 3)
  
  fup_A_4 <- generate_PSA_FUP(baseline_mean_A, baseline_A,
                              response_prop_A, correlation_sd, 4)
  
  fup_A_6 <- generate_PSA_FUP(baseline_mean_A, baseline_A,
                              response_prop_A, correlation_sd, 6)
  
  fup_A_8 <- generate_PSA_FUP(baseline_mean_A, baseline_A,
                              response_prop_A, correlation_sd, 8)
    
  fup_A_9 <- generate_PSA_FUP(baseline_mean_A, baseline_A,
                              response_prop_A, correlation_sd, 9)
  
  fup_A_12 <- generate_PSA_FUP(baseline_mean_A, baseline_A,
                              response_prop_A, correlation_sd, 12)
  
  # generate grp B follow-up PSA
  fup_B_2 <- generate_PSA_FUP(baseline_mean_B, baseline_B,
                              response_prop_B, correlation_sd, 2)
  
  fup_B_3 <- generate_PSA_FUP(baseline_mean_B, baseline_B,
                            response_prop_B, correlation_sd, 3)
  
  fup_B_4 <- generate_PSA_FUP(baseline_mean_B, baseline_B,
                              response_prop_B, correlation_sd, 4)
  
  fup_B_6 <- generate_PSA_FUP(baseline_mean_B, baseline_B,
                              response_prop_B, correlation_sd, 6)
  
  fup_B_8 <- generate_PSA_FUP(baseline_mean_B, baseline_B,
                              response_prop_B, correlation_sd, 8)
  
  fup_B_9 <- generate_PSA_FUP(baseline_mean_B, baseline_B,
                              response_prop_B, correlation_sd, 9)
  
  fup_B_12 <- generate_PSA_FUP(baseline_mean_B, baseline_B,
                               response_prop_B, correlation_sd, 12)
  
  df_psa_sim <- data.frame(group = groups,
                           id = patient_id,
                           time = time,
                           PSA = c(baseline_A, baseline_B,
                                   fup_A_2, fup_B_2,
                                   fup_A_3, fup_B_3, 
                                   fup_A_4, fup_B_4,
                                   fup_A_6, fup_B_6, 
                                   fup_A_8, fup_B_8,
                                   fup_A_9, fup_B_9,
                                   fup_A_12, fup_B_12))

  return(df_psa_sim)
}

# Generates simulated PSA follow-up values - helper function
generate_PSA_FUP_adj <- function(
    baseline_mean,
    baseline_sim,
    response_prop, 
    correlation_sd,
    month
) {
  psa_FUP <- (0.10 * baseline_sim) + 
    (0.20 * rnorm(30, mean = 0, sd = correlation_sd))
  
  if(baseline_mean == 4) {
    # Generate x number of FUP responses from short curve
    psa_short <- rep(0.90 * (
      predict(mod_low_short_adj, 
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[1])
    
    # Generate from non-response curve
    psa_non <- rep(0.90 * (
      predict(mod_low_non,
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[3])
    
    if (response_prop[2] != 0) {
      # Generate from slow response curve
      psa_slow <- rep(0.90 * (
        predict(mod_low_slow,
                data.frame(x = month, x_inv = 1/(month+ 1)))[[1]]
      ), response_prop[2])
    } else { psa_slow <- c() }
  }
  
  if(baseline_mean == 20) {
    # Generate x number of FUP responses from short curve
    psa_short <- rep(0.90 * (
      predict(mod_med_short_adj, 
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[1])
    
    # Generate from non-response curve
    psa_non <- rep(0.90 * (
      predict(mod_med_non,
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[3])
    
    if (response_prop[2] != 0) {
      # Generate from slow response curve
      psa_slow <- rep(0.90 * (
        predict(mod_med_slow_adj,
                data.frame(x = month, x_inv = 1/(month+ 1)))[[1]]
      ), response_prop[2])
    } else { psa_slow <- c() }
  }
  
  if(baseline_mean == 100) {
    # Generate x number of FUP responses from short curve
    psa_short <- rep(0.90 * (
      predict(mod_high_short_adj, 
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[1])
    
    # Generate from non-response curve
    psa_non <- rep(0.90 * (
      predict(mod_high_non,
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[3])
    
    if (response_prop[2] != 0) {
      # Generate from slow response curve
      psa_slow <- rep(0.90 * (
        predict(mod_high_slow_adj,
                data.frame(x = month, x_inv = 1/(month+ 1)))[[1]]
      ), response_prop[2])
    } else { psa_slow <- c() }
  }
  
  psa_change <- c(psa_short, psa_slow, psa_non) |> sample()
  
  # floor at 0.1 ng/mL, the functional sensitivity of a standard clinical PSA
  # assay: keeps PSA positive so the log-scale tests stay defined. Non-binding
  # at mCRPC PSA levels; ultrasensitive assays are a post-prostatectomy tool.
  return(pmax(psa_FUP + psa_change, 0.1))
}

generate_PSA_FUP <- function(
    baseline_mean,
    baseline_sim,
    response_prop, 
    correlation_sd,
    month
) {
  psa_FUP <- (0.10 * baseline_sim) + 
    (0.20 * rnorm(30, mean = 0, sd = correlation_sd))
  
  if(baseline_mean == 4) {
    # Generate x number of FUP responses from short curve
    psa_short <- rep(0.90 * (
      predict(mod_low_short, 
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[1])
    
    # Generate from non-response curve
    psa_non <- rep(0.90 * (
      predict(mod_low_non,
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[3])
    
    if (response_prop[2] != 0) {
      # Generate from slow response curve
      psa_slow <- rep(0.90 * (
        predict(mod_low_slow,
                data.frame(x = month, x_inv = 1/(month+ 1)))[[1]]
      ), response_prop[2])
    } else { psa_slow <- c() }
  }
  
  if(baseline_mean == 20) {
    # Generate x number of FUP responses from short curve
    psa_short <- rep(0.90 * (
      predict(mod_med_short, 
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[1])
    
    # Generate from non-response curve
    psa_non <- rep(0.90 * (
      predict(mod_med_non,
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[3])
    
    if (response_prop[2] != 0) {
      # Generate from slow response curve
      psa_slow <- rep(0.90 * (
        predict(mod_med_slow,
                data.frame(x = month, x_inv = 1/(month+ 1)))[[1]]
      ), response_prop[2])
    } else { psa_slow <- c() }
  }
  
  if(baseline_mean == 100) {
    # Generate x number of FUP responses from short curve
    psa_short <- rep(0.90 * (
      predict(mod_high_short, 
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[1])
    
    # Generate from non-response curve
    psa_non <- rep(0.90 * (
      predict(mod_high_non,
              data.frame(x = month, x_inv = 1/(month + 1)))[[1]]
    ), response_prop[3])
    
    if (response_prop[2] != 0) {
      # Generate from slow response curve
      psa_slow <- rep(0.90 * (
        predict(mod_high_slow,
                data.frame(x = month, x_inv = 1/(month+ 1)))[[1]]
      ), response_prop[2])
    } else { psa_slow <- c() }
  }
  
  psa_change <- c(psa_short, psa_slow, psa_non)
  
  # floor at 0.1 ng/mL, the functional sensitivity of a standard clinical PSA
  # assay: keeps PSA positive so the log-scale tests stay defined. Non-binding
  # at mCRPC PSA levels; ultrasensitive assays are a post-prostatectomy tool.
  return(pmax(psa_FUP + psa_change, 0.1))
}

fisher_sim_test <- function(df_psa_wide) {
  fish_test <- df_psa_wide |>
    # Select distinct patients
    select(group, id, psa50_achieved) |>
    distinct() |>
    summarise(.by = group,
              psa50_no = sum(!psa50_achieved),
              psa50_yes = sum(psa50_achieved)) |>
    select(psa50_no, psa50_yes) |>
    # Fisher's exact test only works on contingency tables
    as.matrix() |>
    fisher.test()
  
  return(fish_test$p.value)
}


# log-ratio t-test on the log (sympercent) scale. This is still an UNADJUSTED
# change score, so it stays regression-to-the-mean prone; it is kept as a
# comparator to show it underperforms the baseline-adjusted log_ancova_test.
log_percent_test <- function(df_psa_wide) {
  log_change_A <- df_psa_wide |> filter(group == "A") |>
    mutate(log_change = (log(FUP) - log(baseline)) * 100) |>
    pull(log_change)

  log_change_B <- df_psa_wide |> filter(group == "B") |>
    mutate(log_change = (log(FUP) - log(baseline)) * 100) |>
    pull(log_change)

  # perform t.test(var.equal = FALSE) on them
  log_test <- t.test(log_change_A, log_change_B, var.equal = FALSE)

  # return p-value
  return(log_test$p.value)
}

raw_change_test <- function(df_psa_wide) {
  change_A <- df_psa_wide |> filter(group == "A") |>
    select(change)
  
  change_B <- df_psa_wide |> filter(group == "B") |>
    select(change)
  
  # perform t.test(var.equal = FALSE) on them
  test <- t.test(change_A, change_B, var.equal = FALSE)
  
  # return p-value
  return(test$p.value)
}

# perchange_sim_test
perchange_sim_test <- function(df_psa_wide) {
  perchange_A <- df_psa_wide |> filter(group == "A") |>
    mutate(perchange = ((FUP - baseline) / baseline) * 100) |>
    select(perchange)
  
  perchange_B <- df_psa_wide |> filter(group == "B") |>
    mutate(perchange = ((FUP - baseline) / baseline) * 100) |>
    select(perchange)
  
  test <- t.test(perchange_A, perchange_B, var.equal = FALSE)
  
  return(test$p.value)
}


# ancova_test
ancova_test <- function(df_psa_wide) {
  mod <- aov(FUP ~ baseline + group, data = df_psa_wide)

  # return p-value for group effect after controlling for baseline
  return(anova(mod)[[5]][2])
}


# log-scale ANCOVA: the recommended inference method (Vickers 2001,
# Kaiser 1989, Cole 2000). Adjusts follow-up for baseline on the log scale,
# where PSA is closer to normal and the baseline relationship closer to
# linear; the group effect back-transforms to a sympercent (symmetric
# percent change) for interpretation.
log_ancova_test <- function(df_psa_wide) {
  mod <- aov(log(FUP) ~ log(baseline) + group, data = df_psa_wide)

  # return p-value for group effect after adjusting for log baseline
  return(anova(mod)[[5]][2])
}


# calculate raw AUC per patient then averages
raw_auc_calc <- function(df_psa_sim) {
  # establish sum of AUCs for later averaging
  auc_total <- 0
  
  # go by patient id
  for (i in 1:60) {
    df_patient <- df_psa_sim |> 
      filter(id == i)
    
    auc_current <- fda.usc::int.simpson2(df_patient$time,
                                df_patient$PSA,
                                equi = FALSE,
                                method = "ESR")
    
    auc_total <- auc_total + auc_current
  }
  
  return(auc_total / 60)
}

# calculate normalized AUC per patient then averages
norm_auc_calc <- function(df_psa_sim) {
  # establish sum of AUCs for later averaging
  auc_total <- 0
  
  # go by patient id
  for (i in 1:60) {
    df_patient <- df_psa_sim |> 
      filter(id == i)
    
    auc_current <- fda.usc::int.simpson2(df_patient$time,
                                df_patient$PSA,
                                equi = FALSE,
                                method = "ESR")
    
    auc_max <- df_patient |> filter(time == 0) |> pull(PSA) * 12
    
    norm_auc_current <- auc_current / auc_max
    
    auc_total <- auc_total + norm_auc_current
  }
  
  return(auc_total / 60)
}


get_fit_alpha <- function(df_patient) {
  # try to fit exponential decay to PSA values
  # if too little values to fit - return NA
  tryCatch(
    {
      fit <- nls(PSA ~ SSasymp(time, yf, y0, log_alpha), 
                 data = df_patient)
      return(log(2) / exp(coef(fit)[3]))
    },
    error = function(e){
      return(NA)
    }
  )
}

halflife_calc <- function(df_psa_sim) {
  halflife <- 0
  num_patients <- 0
  
  for (i in 1:60) {
    df_patient <- df_psa_sim |>
      filter(id == i) |>
      mutate(log_conc = log(PSA))
    
    # subset to only the terminal decreasing part of curve
    df_patient <- df_patient[1:which(df_patient$log_conc ==
                                       min(df_patient$log_conc)),]
    
    halflife_curr <- get_fit_alpha(df_patient)
    if (!is.na(halflife_curr)) {
      halflife <- halflife + halflife_curr
      num_patients <- num_patients + 1
    }
  }
  
  return(halflife / num_patients)
}
