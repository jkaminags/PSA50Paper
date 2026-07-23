# Produces data/sim_results_robustness.RDS. Run from simulations/:
#   Rscript run_robustness.R
#
# Adversarial robustness for the efficiency result. The concern is that a clean
# log-normal model with a multiplicative effect makes log-scale ANCOVA optimal
# by construction. We perturb toward the skeptic's case along two axes:
#   effect scale: multiplicative (log-linear) vs additive on the raw ng/mL scale
#   error tails:  Gaussian vs heavy-tailed (t_3)
# and report size-adjusted power for all six tests at a fixed sample size. The
# invariant to check is whether log-scale ANCOVA stays most efficient, and stays
# ahead of PSA50, once the model no longer favors it by assumption.

suppressMessages(library(tidyverse))

LOG_HALF <- log(0.5)
TESTS <- c("fisher", "perchange", "rawchange", "ancova", "log_ratio", "log_ancova")
alpha <- 0.10
M     <- 4000
n     <- 110
mu_b  <- log(20); s_b <- 1.0; decline <- LOG_HALF; s_f <- 1.0; rho <- 0.6

# One replication under a chosen perturbation.
#   effect "mult": treatment shifts mean log follow-up by mult_delta
#   effect "add" : treatment shifts raw follow-up by raw_shift ng/mL
#   err "t": follow-up error is heavy-tailed (t_nu, scaled to the same sd)
sim_once_r <- function(effect, mult_delta, raw_shift, err = "norm", nu = 3) {
  mu_f <- mu_b + decline
  eps_sd <- s_f * sqrt(1 - rho^2)
  draw_eps <- function(nn) if (err == "norm") rnorm(nn, 0, eps_sd)
                           else rt(nn, nu) * (eps_sd / sqrt(nu / (nu - 2)))
  gen_arm <- function(sm, sr) {
    lb <- rnorm(n, mu_b, s_b)
    lf <- mu_f + rho * (s_f / s_b) * (lb - mu_b) + sm + draw_eps(n)
    tibble(baseline = pmax(exp(lb), 0.1), FUP = pmax(exp(lf) + sr, 0.1))
  }
  A <- gen_arm(0, 0)
  B <- if (effect == "mult") gen_arm(mult_delta, 0) else gen_arm(0, raw_shift)
  grp <- factor(rep(c("A", "B"), each = n))
  baseline <- c(A$baseline, B$baseline); FUP <- c(A$FUP, B$FUP)
  change <- FUP - baseline; perchange <- change / baseline * 100
  logratio <- (log(FUP) - log(baseline)) * 100
  psa50 <- FUP < 0.5 * baseline; isA <- grp == "A"
  tab <- matrix(c(sum(psa50[isA]), sum(!psa50[isA]),
                  sum(psa50[!isA]), sum(!psa50[!isA])), 2, byrow = TRUE)
  c(fisher = fisher.test(tab)$p.value,
    perchange = t.test(perchange[isA], perchange[!isA], var.equal = FALSE)$p.value,
    rawchange = t.test(change[isA], change[!isA], var.equal = FALSE)$p.value,
    ancova = anova(lm(FUP ~ baseline + grp))$`Pr(>F)`[2],
    log_ratio = t.test(logratio[isA], logratio[!isA], var.equal = FALSE)$p.value,
    log_ancova = anova(lm(log(FUP) ~ log(baseline) + grp))$`Pr(>F)`[2])
}

# scenarios: label, effect, magnitude, error. Null shares each scenario's error
# distribution with no treatment shift.
scen <- tribble(
  ~scenario,                  ~effect, ~md,   ~rs, ~err,
  "Multiplicative, Gaussian", "mult",  -0.3,   0,  "norm",
  "Multiplicative, heavy-tailed", "mult", -0.3, 0, "t",
  "Raw-additive, Gaussian",   "add",    0,   -2,  "norm",
  "Raw-additive, heavy-tailed", "add",  0,   -2,  "t"
)

set.seed(20260722)
robustness <- lapply(seq_len(nrow(scen)), function(k) {
  s <- scen[k, ]
  null_p <- t(replicate(M, sim_once_r(s$effect, 0, 0, s$err)))
  alt_p  <- t(replicate(M, sim_once_r(s$effect, s$md, s$rs, s$err)))
  thr <- apply(null_p, 2, quantile, probs = alpha, type = 1, names = FALSE)
  tibble(scenario = s$scenario, test = TESTS,
         cal_size = colMeans(sweep(null_p, 2, thr, "<="))[TESTS],
         cal_power = colMeans(sweep(alt_p, 2, thr, "<="))[TESTS])
}) |> bind_rows() |>
  mutate(scenario = factor(scenario, levels = scen$scenario))

saveRDS(list(robustness = robustness, alpha = alpha, M = M, n = n),
        "data/sim_results_robustness.RDS")

print(robustness |> select(scenario, test, cal_power) |>
        pivot_wider(names_from = test, values_from = cal_power) |>
        as.data.frame(), digits = 3)
