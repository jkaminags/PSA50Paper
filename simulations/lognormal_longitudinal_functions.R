# Log-linear random-slope engine for the PSA doubling time (progression) comparison.
#
# The response engine (lognormal_functions.R) reads a level contrast off a two-
# timepoint trajectory; this one reads a slope off a multi-timepoint trajectory.
# Both are functionals of the same primitive, the log-PSA curve:
#
#   log(PSA_ij) = a_i + b_i * t_ij + e_ij
#
# a_i is the patient's log baseline, b_i the log-scale growth rate, and the
# doubling time is PSADT_i = log(2) / b_i. Treatment slows the rise by shifting
# mean b, so the estimand is a growth-rate difference that back-transforms to a
# doubling-time ratio: the progression-side analog of the sympercent.
#
# The setting is rising PSA (biochemical recurrence / nmCRPC), not the declining
# PSA of the response setting. Slopes are positive and treatment makes them
# smaller.
#
# Note the mechanism differs from the response side. Percentage change loses
# efficiency by dividing by a noisy baseline (regression to the mean); PSADT
# loses it to the nonlinearity of the reciprocal and to per-patient slopes
# estimated with unequal precision. Same pathology, different proximate cause.

suppressMessages(library(tidyverse))

# Non-rising patients have a negative or undefined doubling time. Practice caps
# them at a long ceiling ("effectively not doubling") and we do the same, at 60
# months. The cap compresses the upper tail, which lowers the variance of the
# PSADT t-test, so the convention is conservative toward the endpoint we are
# criticizing. That the rule is needed at all is part of the case against the
# reciprocal scale.
DT_CAP <- 60

dt_from_slope <- function(b) ifelse(b <= log(2) / DT_CAP, DT_CAP, log(2) / b)

# Doubling time implied by a mean growth rate, and the growth rate implied by a
# doubling time. Used to set the treatment effect in clinical units.
dt_of  <- function(b)  log(2) / b
rate_of <- function(dt) log(2) / dt

TESTS_LONG <- c("fisher_dt", "psadt_t", "log_psadt_t", "raw_slope_t", "slope_t")

# One replication on a balanced visit schedule: p-values for the competing
# progression analyses, n patients per arm.
#   delta_b  treatment shift in mean growth rate (negative = slower rise)
#   times    visit times in months, shared by all patients
#   mu_a,s_a mean / sd of log baseline PSA
#   b0       control-arm mean growth rate per month
#   s_b      between-patient sd of the growth rate
#   s_e      within-patient (assay + biological) sd on the log scale
#   dt_cut   rapid-vs-slow doubling-time threshold in months
sim_once_long <- function(n, delta_b, times, mu_a, s_a, b0, s_b, s_e, dt_cut,
                          err = "norm", nu = 3) {
  J <- length(times)
  tc <- times - mean(times)
  Sxx <- sum(tc^2)

  gen_arm <- function(nn, shift) {
    a <- rnorm(nn, mu_a, s_a)
    b <- rnorm(nn, b0 + shift, s_b)
    E <- if (err == "norm") matrix(rnorm(nn * J, 0, s_e), nn, J)
         else matrix(rt(nn * J, nu) * (s_e / sqrt(nu / (nu - 2))), nn, J)
    # floor at 0.1 ng/mL (assay functional sensitivity)
    pmax(outer(a, rep(1, J)) + outer(b, times) + E, log(0.1))
  }

  logY <- rbind(gen_arm(n, 0), gen_arm(n, delta_b))
  isA  <- rep(c(TRUE, FALSE), each = n)

  # balanced schedule, so every patient's OLS slope is a common linear contrast
  slope     <- as.vector(logY %*% tc) / Sxx        # log-scale growth rate
  raw_slope <- as.vector(exp(logY) %*% tc) / Sxx   # raw ng/mL per month
  psadt     <- dt_from_slope(slope)
  rapid     <- psadt < dt_cut

  tab <- matrix(c(sum(rapid[isA]),  sum(!rapid[isA]),
                  sum(rapid[!isA]), sum(!rapid[!isA])),
                nrow = 2, byrow = TRUE)

  tt <- function(x) t.test(x[isA], x[!isA], var.equal = FALSE)$p.value

  c(fisher_dt   = fisher.test(tab)$p.value,
    psadt_t     = tt(psadt),
    log_psadt_t = tt(log(psadt)),
    raw_slope_t = tt(raw_slope),
    slope_t     = tt(slope))
}

# Size-adjusted (calibrated) power at one setting. Each test's rejection
# threshold is its own alpha-quantile under the null, so every test controls
# Type I at exactly alpha and the comparison is not confounded by tests that are
# conservative (Fisher) or liberal at the nominal cutoff. Mirrors cal_power_at()
# in lognormal_functions.R.
cal_power_at_long <- function(n, delta_b, times, mu_a, s_a, b0, s_b, s_e, dt_cut,
                              alpha, M, err = "norm") {
  rep_p <- function(d) t(replicate(M, sim_once_long(n, d, times, mu_a, s_a, b0,
                                                    s_b, s_e, dt_cut, err)))
  null_p <- rep_p(0)
  alt_p  <- rep_p(delta_b)
  thr <- apply(null_p, 2, quantile, probs = alpha, type = 1, names = FALSE)
  list(cal_power = colMeans(sweep(alt_p,  2, thr, "<=")),
       cal_size  = colMeans(sweep(null_p, 2, thr, "<=")),
       raw_power = colMeans(alt_p <= alpha))
}

# ---- unbalanced schedules ----------------------------------------------------
# Under a balanced schedule every per-patient slope has the same variance, so the
# two-stage t-test on slopes is Gauss-Markov optimal and the mixed model can only
# tie it. The mixed model earns its keep when visit counts vary: per-patient
# slopes then carry unequal precision, which the two-stage t-test ignores and the
# mixed model weights correctly. These functions generate dropout and fit both.

# Long-format panel with dropout: patient i contributes visits 1..J_i.
gen_panel <- function(nn, shift, id0, times, mu_a, s_a, b0, s_b, s_e, min_visits) {
  J <- length(times)
  a <- rnorm(nn, mu_a, s_a)
  b <- rnorm(nn, b0 + shift, s_b)
  logY <- pmax(outer(a, rep(1, J)) + outer(b, times) +
                 matrix(rnorm(nn * J, 0, s_e), nn, J), log(0.1))
  # index into the vector rather than sample() it: sample(k:k, ...) would be read
  # as sample(1:k, ...) and silently produce patients with too few visits
  opts <- min_visits:J
  Ji <- opts[sample.int(length(opts), nn, replace = TRUE)]
  tibble(id     = rep(id0 + seq_len(nn), times = J),   # as.vector is column-major
         t      = rep(times, each = nn),
         logpsa = as.vector(logY),
         keep   = as.vector(outer(Ji, seq_len(J), ">="))) |>
    filter(keep) |>
    select(-keep)
}

# One replication under dropout: two-stage slope t-test vs the mixed model.
# Both statistics are monotone in the evidence, and thresholds are calibrated
# under the null, so the normal approximation for the mixed-model p-value costs
# nothing.
sim_once_unbalanced <- function(n, delta_b, times, mu_a, s_a, b0, s_b, s_e,
                                min_visits = 2) {
  d <- bind_rows(gen_panel(n, 0,       0, times, mu_a, s_a, b0, s_b, s_e, min_visits),
                 gen_panel(n, delta_b, n, times, mu_a, s_a, b0, s_b, s_e, min_visits)) |>
    mutate(grp = factor(if_else(id <= n, "A", "B")))

  per_pt <- d |>
    summarise(.by = c(id, grp), slope = cov(t, logpsa) / var(t))
  isA <- per_pt$grp == "A"

  p_two <- t.test(per_pt$slope[isA], per_pt$slope[!isA], var.equal = FALSE)$p.value

  fit <- suppressMessages(suppressWarnings(
    lme4::lmer(logpsa ~ t * grp + (t | id), data = d,
               control = lme4::lmerControl(calc.derivs = FALSE))))
  tval <- summary(fit)$coefficients["t:grpB", "t value"]
  p_lmm <- 2 * pnorm(-abs(tval))

  c(slope_t = p_two, lmm = p_lmm)
}

cal_power_unbalanced <- function(n, delta_b, times, mu_a, s_a, b0, s_b, s_e,
                                 alpha, M, min_visits = 2) {
  rep_p <- function(d) t(replicate(M, sim_once_unbalanced(n, d, times, mu_a, s_a,
                                                          b0, s_b, s_e, min_visits)))
  null_p <- rep_p(0)
  alt_p  <- rep_p(delta_b)
  thr <- apply(null_p, 2, quantile, probs = alpha, type = 1, names = FALSE)
  list(cal_power = colMeans(sweep(alt_p,  2, thr, "<=")),
       cal_size  = colMeans(sweep(null_p, 2, thr, "<=")))
}
