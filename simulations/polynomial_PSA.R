library(tidyverse)

set.seed(1)

# polynomial fits for low baseline PSA value starting at 4 ng/mL

# short responder decreases to 0.5 ng/mL at 4 months
df_low_short <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(8, 4, 3, 2, 1, 0.70, 0.5, 0.55, 0.63, 
        0.90, 1.20, 1.46, 1.83, 2.15, 2.5, 2.9, 3.3)
) |>
  mutate(x_inv = 1/(x+1))

mod_low_short <- lm(y ~ poly(x, 2) + 
                      poly(x_inv, 3), data = df_low_short)

# slow responder decreases to 0.5 ng/mL at 9 months
df_low_slow <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(5, 4, 3.5, 3, 2.3, 1.8, 1.35, 0.95, 0.7, 
        0.58, 0.52, 0.5, 0.6, 0.9, 1.3, 2, 2.7)
) |>
  mutate(x_inv = 1/(x+1))

mod_low_slow <- lm(y ~ poly(x, 2) + 
                     poly(x_inv, 3), data = df_low_slow)

# non responder decreases to 2.5 ng/mL at 3 months
df_low_non <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14),
  y = c(5, 4, 3.5, 3, 2.65, 2.5, 2.35, 2.3, 2.31, 2.34, 2.41, 
        2.53, 2.65, 2.8, 3.0, 3.2, 3.35, 3.6, 3.8)
) |>
  mutate(x_inv = 1/(x+1))

mod_low_non <- lm(y ~ poly(x, 2) + 
                    poly(x_inv, 3), data = df_low_non)

# Medium PSA baseline start, short responder starting at 20 ng/Ml
# decreases to 5 ng/mL at 4 months
df_med_short <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 2, 3, 3.5, 3.5, 4, 4, 4, 4,
        5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(30, 20, 14, 10, 7.5, 6, 5.4, 5.4, 5, 5, 5, 5, 5.3, 6, 
        7.2, 8.4, 9.8, 11.5, 12.8, 14.5, 16, 17.5)
) |>
  mutate(x_inv = 1/(x+1))

mod_med_short <- lm(y ~ poly(x, 2) + 
                      poly(x_inv, 3), data = df_med_short)

# slow responder decreases to 5 ng/mL at 9 months
df_med_slow <- data.frame(
  x = c(-0.5, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(28, 20, 17.5, 15, 13, 10.5, 9, 7.5, 
        6.2, 5.5, 5, 5.2, 6, 7, 8, 10.5)
) |>
  mutate(x_inv = 1/(x+1))

mod_med_slow <- lm(y ~ poly(x, 2) + 
                     poly(x_inv, 3), data = df_med_slow)

# non responder decreases to 5 ng/mL at 3 months
df_med_non <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3, 3, 3.5, 3.5, 3.5,
        4, 5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14),
  y = c(23, 20, 19, 18, 17, 16.2, 15.5, 15, 15, 15, 15.03, 15.03,
        15.03, 15.1, 15.25, 15.5, 
        15.9, 16.3, 16.8, 17.4, 17.8, 18.33, 18.75, 19.3)
) |>
  mutate(x_inv = 1/(x+1))

mod_med_non <- lm(y ~ poly(x, 2) + 
                    poly(x_inv, 3), data = df_med_non)

# High baseline PSA responder curves, starting at 100 ng/mL
# Short responder decreases to 30 ng/mL at 4 months
df_high_short <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 2, 3, 3.5, 3.5, 4, 4, 4, 4,
        5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(120, 100, 90, 75, 55, 40, 35, 35, 30, 30, 30, 30, 
        31, 33.5, 38, 42.5, 48.5, 54, 60, 65.5, 70, 75)
) |>
  mutate(x_inv = 1/(x+1))


mod_high_short <- lm(y ~ poly(x, 3) + 
                       poly(x_inv, 3), data = df_high_short)

# Slow responder decreases to 30 ng/mL at 9 months
df_high_slow <- data.frame(
  x = c(-0.5, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(110, 100, 85, 76, 67, 58, 50, 42, 
        36, 33, 30, 32, 38, 45, 54, 65)
) |>
  mutate(x_inv = 1/(x+1))

mod_high_slow <- lm(y ~ poly(x, 3) + 
                      poly(x_inv, 3), data = df_high_slow)

# Non responder decreases to 70 ng/mL at 3 months
df_high_non <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3, 3, 3.5, 3.5, 3.5,
        4, 5, 6, 7, 8, 9, 
        10, 11, 12, 13, 14),
  y = c(110, 100, 88, 82, 78, 75, 72, 
        70, 70, 70, 70.5, 70.5, 70.5, 71.1, 72.7, 74.8, 
        76.5, 78.8, 80.5, 82.5, 84.5, 87, 90, 92)
) |>
  mutate(x_inv = 1/(x+1))

mod_high_non <- lm(y ~ poly(x, 2) + 
                     poly(x_inv, 3), data = df_high_non)

df_low_short_adj <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(8, 4, 2, 1, 0.75, 0.61, 0.5, 0.4, 0.45,
        0.50, 0.59, 0.70, 0.83, 0.95, 1.1, 1.21, 1.4)
) |>
  mutate(x_inv = 1/(x+1))

mod_low_short_adj <- lm(y ~ poly(x, 2) +
                      poly(x_inv, 3), data = df_low_short_adj)


df_low_slow_adj <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(5, 4, 3.3, 2.8, 2.1, 1.65, 1.31, 0.95, 0.7,
        0.58, 0.52, 0.5, 0.55, 0.66, 0.85, 1.1, 1.5)
) |>
  mutate(x_inv = 1/(x+1))

mod_low_slow_adj <- lm(y ~ poly(x, 2) +
                     poly(x_inv, 3), data = df_low_slow_adj)


df_med_short_adj <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 2, 3, 3.5, 3.5, 4, 4, 4, 4,
        5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(30, 20, 12.5, 8, 5.5, 4.3, 4, 4, 3.8, 3.8, 3.8, 3.8, 4,
        4.15, 4.3, 4.48, 4.95, 5.35, 5.85, 6.3, 6.95, 7.6)
) |>
  mutate(x_inv = 1/(x+1))

mod_med_short_adj <- lm(y ~ poly(x, 2) +
                      poly(x_inv, 3), data = df_med_short_adj)

df_med_slow_adj <- data.frame(
  x = c(-0.5, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(28, 20, 14.5, 11, 9, 7.2, 6.1, 5,
        4.55, 4.15, 4, 4.1, 4.7, 5.5, 6.6, 8)
) |>
  mutate(x_inv = 1/(x+1))

mod_med_slow_adj <- lm(y ~ poly(x, 2) +
                     poly(x_inv, 3), data = df_med_slow_adj)

df_high_short_adj <- data.frame(
  x = c(-0.5, 0, 0.5, 1, 2, 3, 3.5, 3.5, 4, 4, 4, 4,
        5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(120, 100, 70, 49, 33, 28, 26, 26, 25, 25, 25, 25,
        25.5, 26.3, 27.5, 29, 30.5, 32, 34, 36.5, 38, 40)
) |>
  mutate(x_inv = 1/(x+1))


mod_high_short_adj <- lm(y ~ poly(x, 3) +
                       poly(x_inv, 3), data = df_high_short_adj)

df_high_slow_adj <- data.frame(
  x = c(-0.5, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  y = c(110, 100, 70, 56, 48, 40.5, 34, 30,
        27, 25.5, 25, 26.5, 28.6, 32, 36, 40)
) |>
  mutate(x_inv = 1/(x+1))

mod_high_slow_adj <- lm(y ~ poly(x, 3) +
                      poly(x_inv, 3), data = df_high_slow_adj)


#############################################################
# Code below was run for the purpose of finding appropriate #
# SDs for the noise to induce desired correlations between  #
# baseline PSA and change                                   # 
#############################################################

# mod_low_short correlation SDs
# 0.18 for 0.80 corr, 0.32 for 0.60 corr, 0.55 for 0.40, 1.2 for 0.2
# low_short_corr <- c()
# for (i in 1:1000) {
#   baseline <- rnorm(30, mean = 4, sd = 0.5)
# 
#   FUP <- 0.90 * (
#     predict(mod_low_short_adj, newdata = data.frame(x = 4,
#                                                 x_inv = 1/5))[[1]]
#   ) + (0.10 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 1.2))
#   FUP <- ifelse(FUP < 0, 0, FUP)
# 
#   low_short_corr <- c(low_short_corr, cor(baseline, FUP))
# }
# 
# sum(low_short_corr) / 1000
# 
# # Low slow correlation SDs
# # 0.18 for 0.8, 0.32 for 0.6, 0.55 for 0.4, 1.2 for 0.2
# low_slow_corr <- c()
# for (i in 1:1000) {
#   baseline <- rnorm(30, mean = 4, sd = 0.5)
#   
#   FUP <- 0.90 * (
#     predict(mod_low_slow, newdata = data.frame(x = 9,
#                                                x_inv = 1/10))[[1]]
#   ) + (0.10 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 0.18))
#   FUP <- ifelse(FUP < 0, 0, FUP)
#   
#   low_slow_corr <- c(low_slow_corr, cor(baseline, FUP))
# }
# 
# sum(low_slow_corr) / 1000
# 
# # Low non correlation SDs
# # 0.18 for 0.8, 0.32 for 0.6, 0.55 for 0.4, 1.2 for 0.2
# low_non_corr <- c()
# for (i in 1:1000) {
#   baseline <- rnorm(30, mean = 4, sd = 0.5)
#   
#   FUP <- 0.90 * (
#     predict(mod_low_non, newdata = data.frame(x = 12,
#                                               x_inv = 1/13))[[1]]
#   ) + (0.10 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 1.2))
#   FUP <- ifelse(FUP < 0, 0, FUP)
#   
#   low_non_corr <- c(low_non_corr, cor(baseline, FUP))
# }
# 
# sum(low_non_corr) / 1000
# 
# # Medium short correlation SDs
# # 1.1 for 0.8, 2 for 0.6, 3.3 for 0.40, 7 for 0.20
# med_short_corr <- c()
# for (i in 1:1000) {
#   baseline <- rnorm(30, mean = 20, sd = 3)
#   
#   FUP <- 0.90 * (
#     predict(mod_med_short, newdata = data.frame(x = 4,
#                                                 x_inv = 1/5))[[1]]
#   ) + (0.1 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 7))
#   FUP <- ifelse(FUP < 0, 0, FUP)
#   
#   med_short_corr <- c(med_short_corr, cor(baseline, FUP))
# }
# 
# sum(med_short_corr) / 1000
# 
# # Medium slow correlation SDs
# # 1.1 for 0.8, 2 for 0.6, 3.3 for 0.40, 7 for 0.20
# med_slow_corr <- c()
# for (i in 1:1000) {
#   baseline <- rnorm(30, mean = 20, sd = 3)
#   
#   FUP <- 0.90 * (
#     predict(mod_med_slow, newdata = data.frame(x = 9,
#                                                x_inv = 1/10))[[1]]
#   ) + (0.10 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 7))
#   FUP <- ifelse(FUP < 0, 0, FUP)
#   
#   med_slow_corr <- c(med_slow_corr, cor(baseline, FUP))
# }
# 
# sum(med_slow_corr) / 1000
# 
# # Medium non correlation SDs
# # 1.1 for 0.8, 2 for 0.6, 3.3 for 0.40, 7 for 0.20
# med_non_corr <- c()
# for (i in 1:1000) {
#   baseline <- rnorm(30, mean = 20, sd = 3)
#   
#   FUP <- 0.90 * (
#     predict(mod_med_non, newdata = data.frame(x = 3,
#                                               x_inv = 1/4))[[1]]
#   ) + (0.10 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 7))
#   FUP <- ifelse(FUP < 0, 0, FUP)
#   
#   med_non_corr <- c(med_non_corr, cor(baseline, FUP))
# }
# 
# sum(med_non_corr) / 1000
# 
# # High slow correlation SDs
# # 3.7 for 0.8, 6.5 for 0.6, 11.1 for 0.40, 23.5 for 0.20
# high_slow_corr <- c()
# for (i in 1:1000) {
#   baseline <- rnorm(30, mean = 100, sd = 10)
#   
#   FUP <- 0.90 * (
#     predict(mod_high_slow, newdata = data.frame(x = 9,
#                                                 x_inv = 1/10))[[1]]
#   ) + (0.10 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 23.5))
#   FUP <- ifelse(FUP < 0, 0, FUP)
#   
#   high_slow_corr <- c(high_slow_corr, cor(baseline, FUP))
# }
# 
# sum(high_slow_corr) / 1000
# 
# # High slow correlation SDs
# # 3.7 for 0.8, 6.5 for 0.6, 11.1 for 0.40, 23.5 for 0.20
# high_slow_corr <- c()
# for (i in 1:1000) {
#   baseline <- rnorm(30, mean = 100, sd = 10)
#   
#   FUP <- 0.90 * (
#     predict(mod_high_slow, newdata = data.frame(x = 9,
#                                                 x_inv = 1/10))[[1]]
#   ) + (0.10 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 23.5))
#   FUP <- ifelse(FUP < 0, 0, FUP)
#   
#   high_slow_corr <- c(high_slow_corr, cor(baseline, FUP))
# }
# 
# sum(high_slow_corr) / 1000
# 
# # High non correlation SDs
# # 3.7 for 0.8, 6.5 for 0.6, 11.1 for 0.40, 23.5 for 0.20
# high_non_corr <- c()
# for (i in 1:1000) {
#   baseline <- rnorm(30, mean = 100, sd = 10)
#   
#   FUP <- 0.90 * (
#     predict(mod_high_non, newdata = data.frame(x = 3,
#                                                x_inv = 1/4))[[1]]
#   ) + (0.10 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 11.1))
#   FUP <- ifelse(FUP < 0, 0, FUP)
#   
#   high_non_corr <- c(high_non_corr, cor(baseline, FUP))
# }
# 
# sum(high_non_corr) / 1000




## Redefine polynomials based on real data ??
# time <- rep(seq(0:12), each = 30)
# baseline = rnorm(30, 100, 10)
# FUP <- c()
# for (i in 1:12) {
#   FUP_curr <- 0.90 * (predict(mod_high_short_adj, newdata =
#                                 data.frame(x = i,
#                                            x_inv = 1/(i + 1)))) +
#     (0.10 * baseline) + (0.2 * rnorm(30, mean = 0, sd = 23.5))
# 
#   FUP_curr <- ifelse(FUP_curr < 0, 0, FUP_curr)
#   FUP <- c(FUP, FUP_curr)
# }
# dftemp <- data.frame(time = time,
#                      PSA = c(baseline, FUP))
# 
# plot(PSA ~ time, data = dftemp)