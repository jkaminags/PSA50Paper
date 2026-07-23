# Produces data/sim_results_kinetics.RDS, the artifact the paper reads for the
# PSA doubling time (progression) story. Run from the simulations/ directory:
#   Rscript run_kinetics.R
#
# The progression-side counterpart to run_efficiency.R. Same experimental
# apparatus (size-adjusted power, sample size for 80% power, a decomposition
# waterfall), applied to the slope functional of the log-PSA trajectory instead
# of the level functional:
#   (1) calibrated power over sample size   -> N* for 80% power per analysis
#   (2) headline sample-size multiplier for the rapid/slow flag
#   (3) decomposition: dichotomization vs the reciprocal scale
#   (4) kinetics panel (follow-up duration x measurement noise) supporting the
#       scale critique -- the analog of the response side's RTM panel, though the
#       mechanism is the reciprocal's nonlinearity, not baseline division
#   (5) adversarial robustness: heavy-tailed error x follow-up duration
#   (6) mixed model vs two-stage under a balanced schedule and under dropout

source("lognormal_longitudinal_functions.R")

alpha   <- 0.10
M       <- 4000
target  <- 0.80

times   <- c(0, 3, 6, 9, 12)  # quarterly PSA over a year
mu_a    <- log(5)             # baseline median 5 ng/mL (biochemical recurrence)
s_a     <- 0.8
dt0     <- 9                  # control-arm median doubling time, months
dt1     <- 11                 # treated-arm median doubling time
b0      <- rate_of(dt0)
delta_b <- rate_of(dt1) - rate_of(dt0)
s_b     <- 0.03               # between-patient sd of the growth rate
s_e     <- 0.20               # within-patient sd on the log scale
dt_cut  <- dt0                # rapid/slow threshold at the control median,
                              # so exactly half the control arm is flagged: the
                              # point most favorable to the binary test

# ---- (1) calibrated power curve over N ---------------------------------------
set.seed(1)
Ngrid <- c(30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 170, 200, 240, 280, 340, 400)
curve <- lapply(Ngrid, function(n) {
  r <- cal_power_at_long(n, delta_b, times, mu_a, s_a, b0, s_b, s_e, dt_cut, alpha, M)
  tibble(n = n, test = TESTS_LONG, cal_power = r$cal_power[TESTS_LONG],
         cal_size = r$cal_size[TESTS_LONG], raw_power = r$raw_power[TESTS_LONG])
}) |> bind_rows()

nstar_of <- function(df) {
  df <- df |> arrange(n)
  if (max(df$cal_power) < target) return(NA_real_)
  i <- which(df$cal_power >= target)[1]
  if (i == 1) return(df$n[1])
  x0 <- df$n[i - 1]; x1 <- df$n[i]; y0 <- df$cal_power[i - 1]; y1 <- df$cal_power[i]
  x0 + (target - y0) * (x1 - x0) / (y1 - y0)
}
nstar <- curve |> group_by(test) |> group_modify(~ tibble(n_star = nstar_of(.x))) |> ungroup()
get_ns <- function(t) nstar$n_star[nstar$test == t]
ns_flag <- get_ns("fisher_dt");   ns_dt  <- get_ns("psadt_t")
ns_ldt  <- get_ns("log_psadt_t"); ns_sl  <- get_ns("slope_t")

# ---- (2) headline ------------------------------------------------------------
headline <- list(
  target = target, dt0 = dt0, dt1 = dt1, dt_ratio = dt1 / dt0, delta_b = delta_b,
  n_star_flag = ns_flag, n_star_slope = ns_sl,
  multiplier = ns_flag / ns_sl, discarded_frac = 1 - ns_sl / ns_flag
)

# ---- (3) decomposition (patients per arm) ------------------------------------
# flag -> PSADT t-test isolates dichotomization; PSADT -> growth rate isolates
# the reciprocal scale, split by whether logging the doubling time is enough.
decomp <- tibble(
  step = c("dichotomization", "scale", "log_dt", "reciprocal_residual", "total"),
  label = c("Dichotomization (rapid/slow flag to PSADT t-test)",
            "Reciprocal scale (PSADT to growth rate)",
            "  logging the doubling time (PSADT to log PSADT)",
            "  remaining reciprocal penalty (log PSADT to growth rate)",
            "Total (rapid/slow flag to growth rate)"),
  patients_per_arm = c(ns_flag - ns_dt, ns_dt - ns_sl,
                       ns_dt - ns_ldt, ns_ldt - ns_sl, ns_flag - ns_sl)
)

# ---- (4) kinetics panel: follow-up duration x measurement noise ---------------
# Shorter follow-up and noisier assays make per-patient slopes less precise, push
# more patients toward a non-positive estimated slope, and inflate the reciprocal
# transform. This is the scale mechanism on the progression side.
set.seed(2)
sched <- list(`6`  = c(0, 3, 6),
              `12` = c(0, 3, 6, 9, 12),
              `24` = c(0, 3, 6, 9, 12, 18, 24))
panel_grid <- expand.grid(dur = names(sched), s_e = c(0.15, 0.30),
                          stringsAsFactors = FALSE)
kinetics_panel <- lapply(seq_len(nrow(panel_grid)), function(k) {
  tk <- sched[[panel_grid$dur[k]]]
  r <- cal_power_at_long(100, delta_b, tk, mu_a, s_a, b0, s_b,
                         panel_grid$s_e[k], dt_cut, alpha, M)
  tibble(months = as.integer(panel_grid$dur[k]), n_visits = length(tk),
         s_e = panel_grid$s_e[k], test = TESTS_LONG, cal_power = r$cal_power[TESTS_LONG])
}) |> bind_rows()

# fraction of patients whose estimated PSA slope is non-positive, so their
# doubling time is undefined and must be capped
set.seed(3)
undefined_frac <- lapply(seq_len(nrow(panel_grid)), function(k) {
  tk <- sched[[panel_grid$dur[k]]]
  tc <- tk - mean(tk); Sxx <- sum(tc^2); J <- length(tk)
  a <- rnorm(20000, mu_a, s_a); b <- rnorm(20000, b0, s_b)
  logY <- pmax(outer(a, rep(1, J)) + outer(b, tk) +
                 matrix(rnorm(20000 * J, 0, panel_grid$s_e[k]), 20000, J), log(0.1))
  sl <- as.vector(logY %*% tc) / Sxx
  tibble(months = as.integer(panel_grid$dur[k]), s_e = panel_grid$s_e[k],
         frac_capped = mean(sl <= log(2) / DT_CAP))
}) |> bind_rows()

# ---- (5) adversarial robustness ----------------------------------------------
set.seed(20260723)
rob_grid <- expand.grid(dur = c("6", "12", "24"), err = c("norm", "t"),
                        stringsAsFactors = FALSE)
robustness <- lapply(seq_len(nrow(rob_grid)), function(k) {
  tk <- sched[[rob_grid$dur[k]]]
  r <- cal_power_at_long(100, delta_b, tk, mu_a, s_a, b0, s_b, s_e, dt_cut,
                         alpha, M, err = rob_grid$err[k])
  tibble(scenario = paste0(rob_grid$dur[k], " months, ",
                           if (rob_grid$err[k] == "norm") "Gaussian" else "heavy-tailed"),
         months = as.integer(rob_grid$dur[k]), err = rob_grid$err[k],
         test = TESTS_LONG, cal_power = r$cal_power[TESTS_LONG],
         cal_size = r$cal_size[TESTS_LONG])
}) |> bind_rows()

# ---- (6) mixed model vs two-stage, balanced vs dropout -----------------------
# Under a balanced schedule every per-patient slope has the same variance, so the
# two-stage t-test is Gauss-Markov optimal and the mixed model can only tie it.
# Under dropout the slopes carry unequal precision and only the mixed model
# weights them correctly. M is smaller here because each replication fits lmer.
set.seed(11)
M_lmm <- 1500
lmm_panel <- lapply(c(5, 2), function(mv) {
  r <- cal_power_unbalanced(80, delta_b, times, mu_a, s_a, b0, s_b, s_e,
                            alpha, M_lmm, min_visits = mv)
  tibble(schedule = if (mv == 5) "Balanced (5 visits)" else "Dropout (2-5 visits)",
         test = c("slope_t", "lmm"),
         cal_power = r$cal_power[c("slope_t", "lmm")],
         cal_size  = r$cal_size[c("slope_t", "lmm")])
}) |> bind_rows()

saveRDS(list(
  meta = list(alpha = alpha, M = M, M_lmm = M_lmm, times = times, mu_a = mu_a,
              s_a = s_a, b0 = b0, delta_b = delta_b, s_b = s_b, s_e = s_e,
              dt0 = dt0, dt1 = dt1, dt_cut = dt_cut, dt_cap = DT_CAP,
              target = target),
  curve = curve, nstar = nstar, headline = headline, decomp = decomp,
  kinetics_panel = kinetics_panel, undefined_frac = undefined_frac,
  robustness = robustness, lmm_panel = lmm_panel
), "data/sim_results_kinetics.RDS")

cat(sprintf("HEADLINE: rapid/slow flag %.0f/arm vs growth rate %.0f/arm = %.2fx (%.0f%% discarded)\n",
            ns_flag, ns_sl, headline$multiplier, 100 * headline$discarded_frac))
cat(sprintf("N*: flag %.0f | PSADT %.0f | log PSADT %.0f | growth rate %.0f\n",
            ns_flag, ns_dt, ns_ldt, ns_sl))
print(as.data.frame(decomp[, c("step", "patients_per_arm")]), digits = 3)
print(as.data.frame(lmm_panel), digits = 3)
