# Produces data/sim_results_efficiency.RDS, the artifact the paper reads for the
# efficiency story. Run from the simulations/ directory:
#   Rscript run_efficiency.R
#
# Experiment on the clean log-normal model (lognormal_functions.R):
#   (1) calibrated power over sample size  -> N* for 80% power per test
#   (2) sample-size multiplier + discarded fraction (headline)
#   (3) two-critique decomposition in patients/arm (waterfall)
#   (4) RTM panel (rho x baseline heterogeneity) supporting the scale critique
# The mechanism-based robustness ranking is produced by run_robustness.R.

source("lognormal_functions.R")

alpha    <- 0.10
M        <- 4000
mu_b     <- log(20)   # baseline median 20 ng/mL (mCRPC-plausible)
s_b      <- 1.0       # heterogeneous baseline (sd of log PSA)
decline  <- LOG_HALF  # control arm at a 50% median decline
s_f      <- 1.0
rho_base <- 0.6
delta    <- -0.3      # 30 sympercents (~26% lower geometric-mean PSA on treatment)
target   <- 0.80

# ---- (1) calibrated power curve over N ---------------------------------------
set.seed(1)
Ngrid <- c(40, 50, 60, 70, 80, 90, 100, 110, 120, 140, 160, 180, 200, 230, 260, 300, 340)
curve <- lapply(Ngrid, function(n) {
  r <- cal_power_at(n, delta, mu_b, s_b, decline, s_f, rho_base, alpha, M)
  tibble(n = n, test = TESTS, cal_power = r$cal_power[TESTS],
         cal_size = r$cal_size[TESTS], raw_power = r$raw_power[TESTS])
}) |> bind_rows()

# N* = per-arm sample size reaching target power (linear interp on bracket)
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
ns_fisher <- get_ns("fisher"); ns_per <- get_ns("perchange")
ns_lr <- get_ns("log_ratio"); ns_lanc <- get_ns("log_ancova")

# ---- (2) headline ------------------------------------------------------------
headline <- list(
  target = target, delta = delta, sympercent = 100 * delta,
  n_star_fisher = ns_fisher, n_star_logancova = ns_lanc,
  multiplier = ns_fisher / ns_lanc, discarded_frac = 1 - ns_lanc / ns_fisher
)

# ---- (3) two-critique waterfall (patients per arm) ---------------------------
# PSA50 -> % change isolates dichotomization; % change -> log-ANCOVA isolates
# the scale + baseline adjustment; log-ratio splits the latter.
decomp <- tibble(
  step = c("dichotomization", "scale_and_adjustment",
           "scale_swap", "baseline_adjustment", "total"),
  label = c("Dichotomization (PSA50 to % change)",
            "Scale + baseline adjustment (% change to log-ANCOVA)",
            "  scale swap (% change to log ratio)",
            "  baseline adjustment (log ratio to log-ANCOVA)",
            "Total (PSA50 to log-ANCOVA)"),
  patients_per_arm = c(ns_fisher - ns_per, ns_per - ns_lanc,
                       ns_per - ns_lr, ns_lr - ns_lanc, ns_fisher - ns_lanc)
)

# ---- (4) RTM panel -----------------------------------------------------------
set.seed(2)
panel_grid <- expand.grid(rho = c(0.3, 0.6, 0.85), s_b = c(0.6, 1.0))
rtm_panel <- lapply(seq_len(nrow(panel_grid)), function(k) {
  r <- cal_power_at(110, delta, mu_b, panel_grid$s_b[k], decline, s_f,
                    panel_grid$rho[k], alpha, M)
  tibble(rho = panel_grid$rho[k], s_b = panel_grid$s_b[k],
         test = TESTS, cal_power = r$cal_power[TESTS])
}) |> bind_rows()

saveRDS(list(
  meta = list(alpha = alpha, M = M, mu_b = mu_b, s_b = s_b, decline = decline,
              s_f = s_f, rho = rho_base, delta = delta, target = target),
  curve = curve, nstar = nstar, headline = headline,
  decomp = decomp, rtm_panel = rtm_panel
), "data/sim_results_efficiency.RDS")

cat(sprintf("HEADLINE: PSA50 %.0f/arm vs log-ANCOVA %.0f/arm = %.2fx (%.0f%% discarded)\n",
            ns_fisher, ns_lanc, headline$multiplier, 100 * headline$discarded_frac))
