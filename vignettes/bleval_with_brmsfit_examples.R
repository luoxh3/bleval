################################################################################
## Example: Evaluating Bayesian hierarchical models using bleval_with_brmsfit ##
##          for different random-effects configurations (preliminary checks)  ##
##                                                                            ##
## This script demonstrates 6 scenarios of two-level hierarchical models      ##
## fitted with brms, using bleval to compute information criteria             ##
## (DIC, WAIC, LOOIC) and log fully marginal likelihoods.                     ##
##                                                                            ##
## Scenarios:                                                                 ##
## 1. Gaussian: Random intercept only                                         ##
## 2. Gaussian: Random slope only                                             ##
## 3. Gaussian: Correlated random intercept and slope                         ##
## 4. Gaussian: Three correlated random effects (intercept + two slopes)      ##
## 5. Bernoulli: Correlated random intercept and slope                        ##
## 6. Poisson: Correlated random intercept and slope                          ##
##                                                                            ##
## Authors: Xiaohui Luo, Jieyuan Dong                                         ##
## Last modified: 11/11/2025                                                  ##
################################################################################

suppressPackageStartupMessages({
 library(brms)
 library(posterior)
 library(mvtnorm)
 library(Matrix)
 library(bleval)
})

source("bleval_with_brmsfit.R")

set.seed(1107)

make_priors <- function(K, family = "gaussian") {
  # - Always: Normal(0,5) for fixed effects b.
  # - If K>=1: half-Cauchy (implemented as Cauchy truncated at >0 in brms) on sd.
  # - If K>=2: LKJ(2) correlation prior.
  # - Only for gaussian: Cauchy(0,2) for residual sigma.
  pri <- prior(normal(0, 5), class = "b")
  if (K >= 1) pri <- c(pri, prior(cauchy(0, 2), class = "sd"))
  if (K >= 2) pri <- c(pri, prior(lkj(2), class = "cor"))
  if (family == "gaussian") pri <- c(pri, prior(cauchy(0, 2), class = "sigma"))
  pri
}

# ================================================================================
# Scenario 1: Gaussian
# 1 random intercept; fixed slope for x1
# ================================================================================
N <- 100; T_per <- 20
x1 <- rnorm(N * T_per)
group <- factor(rep(seq_len(N), each = T_per))
# K = 1 (Intercept only)
Sigma1 <- matrix(1^2, nrow = 1, dimnames = list("Intercept", NULL))
u_mat <- mvtnorm::rmvnorm(N, mean = rep(0, 1), sigma = Sigma1)
beta <- c(Intercept = 0, x1 = 0.5)
eta <- numeric(N * T_per)
for (g in seq_len(N)) {
  idx_g <- which(group == g)
  ru_int <- u_mat[g, 1]
  eta[idx_g] <- (beta["Intercept"] + ru_int) + beta["x1"] * x1[idx_g]
}
y <- eta + rnorm(N * T_per, 0, 0.5)
sim_dat <- data.frame(y = y, x1 = x1, group = group)
PRIORS1 <- make_priors(K = 1)
fit1 <- brm(formula = bf(y ~ x1 + (1 | group)), data = sim_dat, prior = PRIORS1,
            iter = 5000, warmup = 2500, chains = 2, cores = 2, seed = 1111, refresh = 0)
data_list1 <- list(Nobs = nrow(sim_dat), subject = as.integer(sim_dat$group), N = length(unique(sim_dat$group)),
                   y = sim_dat$y, x1 = sim_dat$x1)
res1 <- bleval_with_brmsfit(fit = fit1, data = data_list1, Ngrid = 3, parallel = FALSE, n_cores = 1)
print(summary(fit1)); print(res1$ic_res); print(res1$fml_res)
#    p_dic     elpd_dic          dic       p_waic    elpd_waic         waic      p_looic   elpd_looic        looic 
# 4.002320 -1674.215151  3348.430301     3.565515 -1674.034488  3348.068975     3.568346 -1674.037319  3348.074638 
# Bridge sampling estimate of the log marginal likelihood: -1688.39
# Estimate obtained in 5 iteration(s) via method "normal".

# ================================================================================
# Scenario 2: Gaussian
# Fixed intercept; random slope on x1 
# ================================================================================
N <- 50; T_per <- 10
x1 <- rnorm(N * T_per)
group <- factor(rep(seq_len(N), each = T_per))
# K = 1 but term name is x1 (no Intercept)
Sigma2 <- matrix(0.55^2, nrow = 1, dimnames = list("x1", NULL))
u_mat <- mvtnorm::rmvnorm(N, mean = rep(0, 1), sigma = Sigma2)
beta <- c(Intercept = 0.3, x1 = 0.65)
eta <- numeric(N * T_per)
for (g in seq_len(N)) {
  idx_g <- which(group == g)
  ru_x1 <- u_mat[g, 1]
  eta[idx_g] <- beta["Intercept"] + (beta["x1"] + ru_x1) * x1[idx_g]
}
y <- eta + rnorm(N * T_per, 0, 0.28)
sim_dat <- data.frame(y = y, x1 = x1, group = group)
# brms fit
PRIORS2 <- make_priors(K = 1)
fit2 <- brm(formula = bf(y ~ x1 + (0 + x1 | group)), data = sim_dat, prior = PRIORS2,
            iter = 5000, warmup = 2500, chains = 2, cores = 2, seed = 111, refresh = 0)
data_list2 <- list(Nobs = nrow(sim_dat), subject = as.integer(sim_dat$group), N = length(unique(sim_dat$group)),
                   y = sim_dat$y, x1 = sim_dat$x1)
# model evaluation
res2 <- bleval_with_brmsfit(fit = fit2, data = data_list2, Ngrid = 3, parallel = FALSE, n_cores = 1)
print(summary(fit2)); print(res2$ic_res); print(res2$fml_res)
#    p_dic    elpd_dic         dic      p_waic   elpd_waic        waic     p_looic  elpd_looic       looic 
# 3.959516 -166.761914  333.523828    4.112869 -166.918348  333.836696    4.128951 -166.934429  333.868858 
# Bridge sampling estimate of the log marginal likelihood: -180.748
# Estimate obtained in 5 iteration(s) via method "normal".

# ================================================================================
# Scenario 3: Gaussian
# Correlated random effects (K=2): Random intercept + random slope on x1
# ================================================================================
N <- 50; T_per <- 10
x1 <- rnorm(N * T_per)
group <- factor(rep(seq_len(N), each = T_per))
Sigma3 <- matrix(c(0.5^2, 0.15*0.5*0.4, 0.15*0.5*0.4, 0.4^2), nrow = 2, byrow = TRUE,
                 dimnames = list(c("Intercept","x1"), c("Intercept","x1")))
u_mat <- mvtnorm::rmvnorm(N, mean = rep(0, 2), sigma = Sigma3)
beta <- c(Intercept = 0, x1 = 0.3)
eta <- numeric(N * T_per)
for (g in seq_len(N)) {
  idx_g <- which(group == g)
  ru_int <- u_mat[g, 1]; ru_x1 <- u_mat[g, 2]
  eta[idx_g] <- (beta["Intercept"] + ru_int) + (beta["x1"] + ru_x1) * x1[idx_g]
}
y <- eta + rnorm(N * T_per, 0, 0.25)
sim_dat <- data.frame(y = y, x1 = x1, group = group)
# brms fit
PRIORS3 <- make_priors(K = 2)
fit3 <- brm(formula = bf(y ~ x1 + (x1 | group)), data = sim_dat, prior = PRIORS3,
            iter = 5000, warmup = 2500, chains = 2, cores = 2, seed = 111, refresh = 0)
data_list3 <- list(Nobs = nrow(sim_dat), subject = as.integer(sim_dat$group), N = length(unique(sim_dat$group)),
                   y = sim_dat$y, x1 = sim_dat$x1)
# model evaluation
res3 <- bleval_with_brmsfit(fit = fit3, data = data_list3, Ngrid = 3, parallel = FALSE, n_cores = 1)
print(summary(fit3)); print(res3$ic_res); print(res3$fml_res)
#    p_dic    elpd_dic         dic      p_waic   elpd_waic        waic     p_looic  elpd_looic       looic 
# 5.814436 -169.895435  339.790869    6.332505 -170.343579  340.687158    6.398250 -170.409323  340.818647 
# Bridge sampling estimate of the log marginal likelihood: -184.7061
# Estimate obtained in 6 iteration(s) via method "normal".

# ================================================================================
# Scenario 4: Gaussian
# Correlated random effects (K=3): Random intercept + two random slopes (x1, x2) 
# ================================================================================
N <- 50; T_per <- 10
x1 <- rnorm(N * T_per); x2 <- rnorm(N * T_per)
group <- factor(rep(seq_len(N), each = T_per))
# random effects
sd_u0 <- 1.0
sd_u1 <- 0.5
sd_u2 <- 0.4
rho01 <- 0.2
rho02 <- -0.1
rho12 <- 0.3
Sigma4 <- matrix(c(
  sd_u0^2, sd_u0*sd_u1*rho01, sd_u0*sd_u2*rho02,
  sd_u0*sd_u1*rho01, sd_u1^2, sd_u1*sd_u2*rho12,
  sd_u0*sd_u2*rho02, sd_u1*sd_u2*rho12, sd_u2^2
), nrow = 3, byrow = TRUE,
dimnames = list(c("Intercept","x1","x2"), c("Intercept","x1","x2")))
u_mat <- mvtnorm::rmvnorm(N, mean = rep(0, 3), sigma = Sigma4)
# fixed effects
beta <- c(Intercept = 0, x1 = 0.5, x2 = -0.4) 
eta <- numeric(N * T_per)
for (g in seq_len(N)) {
  idx_g <- which(group == g)
  ru_int <- u_mat[g, 1]; ru_x1 <- u_mat[g, 2]; ru_x2 <- u_mat[g, 3]
  eta[idx_g] <- (beta["Intercept"] + ru_int) + (beta["x1"] + ru_x1) * x1[idx_g] + (beta["x2"] + ru_x2) * x2[idx_g]
}
y <- eta + rnorm(N * T_per, 0, 0.12)
sim_dat <- data.frame(y = y, x1 = x1, x2 = x2, group = group)
# brms fit
PRIORS4 <- make_priors(K = 3)
fit4 <- brm(formula = bf(y ~ x1 + x2 + (x1 + x2 | group)), data = sim_dat, prior = PRIORS4,
            iter = 5000, warmup = 2500, chains = 2, cores = 2, seed = 111, refresh = 0)
data_list4 <- list(Nobs = nrow(sim_dat), subject = as.integer(sim_dat$group), N = length(unique(sim_dat$group)),
                   y = sim_dat$y, x1 = sim_dat$x1, x2 = sim_dat$x2)
# model evaluation
res4 <- bleval_with_brmsfit(fit = fit4, data = data_list4, Ngrid = 3, parallel = FALSE, n_cores = 1)
print(summary(fit4)); print(res4$ic_res); print(res4$fml_res)
#    p_dic   elpd_dic        dic     p_waic  elpd_waic       waic    p_looic elpd_looic      looic 
# 9.303342 -22.966773  45.933546   8.193421 -22.672991  45.345982   8.246736 -22.726306  45.452612 
# Bridge sampling estimate of the log marginal likelihood: -43.50208
# Estimate obtained in 6 iteration(s) via method "normal".

# ================================================================================
# Scenario 5: Bernoulli (Logistic Regression)
# Correlated random effects (K=2): Random intercept + random slope on x1
# ================================================================================
N <- 200; T_per <- 25
x1 <- rnorm(N * T_per)
group <- factor(rep(seq_len(N), each = T_per))
sd_u0 <- 0.5
sd_u1 <- 0.4
rho01 <- 0.3
Sigma5 <- matrix(c(
  sd_u0^2, sd_u0*sd_u1*rho01,
  sd_u0*sd_u1*rho01, sd_u1^2
), nrow = 2, byrow = TRUE,
dimnames = list(c("Intercept","x1"), c("Intercept","x1")))
u_mat <- mvtnorm::rmvnorm(N, mean = rep(0, 2), sigma = Sigma5)
beta <- c(Intercept = 0, x1 = 0.6)
eta <- numeric(N * T_per)
for (g in seq_len(N)) {
  idx_g <- which(group == g)
  ru_int <- u_mat[g, 1]; ru_x1 <- u_mat[g, 2]
  eta[idx_g] <- (beta["Intercept"] + ru_int) + (beta["x1"] + ru_x1) * x1[idx_g]
}
prob <- plogis(eta)
y <- rbinom(N * T_per, size = 1, prob = prob)
sim_dat_logistic <- data.frame(y = y, x1 = x1, group = group)
# brms fit
PRIORS5 <- make_priors(K = 2, "bernoulli")
fit5 <- brm(formula = bf(y ~ x1 + (x1 | group), family = bernoulli()),
            data = sim_dat_logistic, prior = PRIORS5,
            iter = 6000, warmup = 4000, chains = 4, cores = 4, seed = 111, refresh = 0)
data_list5 <- list(Nobs = nrow(sim_dat_logistic), subject = as.integer(sim_dat_logistic$group),
                   N = length(unique(sim_dat_logistic$group)),
                   y = sim_dat_logistic$y, x1 = sim_dat_logistic$x1)
# model evaluation
res5 <- bleval_with_brmsfit(fit = fit5, data = data_list5, Ngrid = 9, parallel = FALSE, n_cores = 1)
print(summary(fit5)); print(res5$ic_res); print(res5$fml_res)
#    p_dic     elpd_dic          dic       p_waic    elpd_waic         waic      p_looic   elpd_looic        looic 
# 4.917250 -3231.519953  6463.039905     4.788961 -3231.488273  6462.976546     4.794551 -3231.493863  6462.987726 
# Bridge sampling estimate of the log marginal likelihood: -3242.39
# Estimate obtained in 5 iteration(s) via method "normal".

# ================================================================================
# Scenario 6: Poisson
# Correlated random effects (K=2): Random intercept + random slope on x1
# ================================================================================
N <- 100; T_per <- 20
x1 <- rnorm(N * T_per)
group <- factor(rep(seq_len(N), each = T_per))
sd_u0 <- 0.5
sd_u1 <- 0.3
rho01 <- 0.2
Sigma6 <- matrix(c(
  sd_u0^2, sd_u0*sd_u1*rho01,
  sd_u0*sd_u1*rho01, sd_u1^2
), nrow = 2, byrow = TRUE,
dimnames = list(c("Intercept","x1"), c("Intercept","x1")))
u_mat <- mvtnorm::rmvnorm(N, mean = rep(0, 2), sigma = Sigma6)
beta <- c(Intercept = 0.3, x1 = 0.5)
eta <- numeric(N * T_per)
for (g in seq_len(N)) {
  idx_g <- which(group == g)
  ru_int <- u_mat[g, 1]; ru_x1 <- u_mat[g, 2]
  eta[idx_g] <- (beta["Intercept"] + ru_int) + (beta["x1"] + ru_x1) * x1[idx_g]
}
lambda <- exp(eta)
y <- rpois(N * T_per, lambda = lambda)
sim_dat_poisson <- data.frame(y = y, x1 = x1, group = group)
# brms fit
PRIORS6 <- make_priors(K = 2, family = "poisson")
fit6 <- brm(formula = bf(y ~ x1 + (x1 | group), family = poisson()), 
            data = sim_dat_poisson, prior = PRIORS6,
            iter = 5000, warmup = 2500, chains = 2, cores = 2, seed = 111, refresh = 0)
data_list6 <- list(Nobs = nrow(sim_dat_poisson), subject = as.integer(sim_dat_poisson$group), 
                   N = length(unique(sim_dat_poisson$group)),
                   y = sim_dat_poisson$y, x1 = sim_dat_poisson$x1)
# model evaluation
res6 <- bleval_with_brmsfit(fit = fit6, data = data_list6, Ngrid = 5, parallel = FALSE, n_cores = 1)
print(summary(fit6)); print(res6$ic_res); print(res6$fml_res)
#    p_dic     elpd_dic          dic       p_waic    elpd_waic         waic      p_looic   elpd_looic        looic 
# 4.854908 -3173.995993  6347.991985     5.604806 -3174.534050  6349.068100     5.660886 -3174.590130  6349.180261 
# Bridge sampling estimate of the log marginal likelihood: -3186.243
# Estimate obtained in 5 iteration(s) via method "normal".
