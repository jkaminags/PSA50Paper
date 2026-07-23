# Clean log-normal generative engine for the PSA50 efficiency comparison.
#
# PSA is log-normal, so (log baseline, log follow-up) is generated as bivariate
# normal with an explicit baseline-follow-up correlation rho. Treatment shifts
# the mean of log follow-up by delta, a log-scale effect whose 100x value is a
# sympercent (Cole 2000). This is the standard trial-design model (Vickers
# 2001); log-scale ANCOVA is efficient under it by construction, which is why
# the mechanism-based response-curve model (sim_functions.R) is run separately
# as a robustness anchor that favors no test.
#
# Baseline median is exp(mu_b); the control arm sits at a chosen median decline
# so the PSA50 flag is informative (neither everyone nor no one responds).

suppressMessages(library(tidyverse))

LOG_HALF <- log(0.5)  # a 50% median decline in the control arm

# One replication: the six competing test p-values for a two-arm trial of n/arm.
#   delta   log-scale treatment effect (negative = deeper PSA decline)
#   mu_b    mean log baseline;   s_b  sd of log baseline
#   decline mean log(follow-up / baseline) in the control arm
#   s_f     sd of log follow-up; rho  baseline-follow-up correlation
sim_once <- function(n, delta, mu_b, s_b, decline, s_f, rho) {
  mu_f <- mu_b + decline
  eps_sd <- s_f * sqrt(1 - rho^2)

  gen_arm <- function(nn, shift) {
    lb <- rnorm(nn, mu_b, s_b)
    lf <- mu_f + rho * (s_f / s_b) * (lb - mu_b) + shift + rnorm(nn, 0, eps_sd)
    # floor at 0.1 ng/mL (assay functional sensitivity); non-binding here
    tibble(baseline = pmax(exp(lb), 0.1), FUP = pmax(exp(lf), 0.1))
  }

  A <- gen_arm(n, 0)
  B <- gen_arm(n, delta)
  grp      <- factor(rep(c("A", "B"), each = n))
  baseline <- c(A$baseline, B$baseline)
  FUP      <- c(A$FUP, B$FUP)

  change    <- FUP - baseline
  perchange <- change / baseline * 100
  logratio  <- (log(FUP) - log(baseline)) * 100
  psa50     <- FUP < 0.5 * baseline
  isA       <- grp == "A"

  tab <- matrix(c(sum(psa50[isA]),  sum(!psa50[isA]),
                  sum(psa50[!isA]), sum(!psa50[!isA])),
                nrow = 2, byrow = TRUE)

  c(
    fisher     = fisher.test(tab)$p.value,
    perchange  = t.test(perchange[isA], perchange[!isA], var.equal = FALSE)$p.value,
    rawchange  = t.test(change[isA],    change[!isA],    var.equal = FALSE)$p.value,
    # group entered last, so its sequential test is baseline-adjusted (ANCOVA)
    ancova     = anova(lm(FUP ~ baseline + grp))$`Pr(>F)`[2],
    log_ratio  = t.test(logratio[isA],  logratio[!isA],  var.equal = FALSE)$p.value,
    log_ancova = anova(lm(log(FUP) ~ log(baseline) + grp))$`Pr(>F)`[2]
  )
}

TESTS <- c("fisher", "perchange", "rawchange", "ancova", "log_ratio", "log_ancova")

# Size-adjusted (calibrated) power at one setting. Each test's rejection
# threshold is its own alpha-quantile under the null, so every test controls
# Type I at exactly alpha and the power comparison is not confounded by tests
# that happen to be conservative (Fisher) or liberal at the nominal cutoff.
cal_power_at <- function(n, delta, mu_b, s_b, decline, s_f, rho, alpha, M) {
  null_p <- t(replicate(M, sim_once(n, 0,     mu_b, s_b, decline, s_f, rho)))
  alt_p  <- t(replicate(M, sim_once(n, delta, mu_b, s_b, decline, s_f, rho)))
  thr    <- apply(null_p, 2, quantile, probs = alpha, type = 1, names = FALSE)
  list(cal_power = colMeans(sweep(alt_p, 2, thr, "<=")),
       cal_size  = colMeans(sweep(null_p, 2, thr, "<=")),
       raw_power = colMeans(alt_p <= alpha),
       raw_size  = colMeans(null_p <= alpha))
}
