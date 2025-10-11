## This R script demonstrates Bayesian evaluation of a location-scale model with
## random effects for both mean and variance structures using the bleval package.
##
## Key characteristics:
## - Model category: Multilevel Model (MLM)
## - Number of latent variables: 4
## - Random intercept and slope for location (mean structure)
## - Random intercept and slope for scale (variance structure)
## - Heteroscedastic errors where residual variance depends on predictors
## - Log-linear link for variance components to ensure positivity

## 1 Data Generation ###########################################################

# Load required packages
library(MASS)
library(dplyr)

set.seed(123)

# Model parameters
Nnum <- 200  # Number of level-2 units (e.g., individuals)
Tnum <- 50   # Number of level-1 units (e.g., time points)

# Location parameters (mean structure)
location_beta <- c(0, 0.3)  # Fixed effects for location: intercept and slope
location_rho <- 0.3         # Correlation between location random intercept and slope
location_u_cov <- sqrt(1) * sqrt(0.09) * location_rho
location_u_sigma <- matrix(c(1, location_u_cov,
                             location_u_cov, 0.09), nrow = 2)

# Scale parameters (variance structure)
scale_beta <- c(0, 0.2)     # Fixed effects for scale: intercept and slope for log(variance)
scale_rho <- 0.2            # Correlation between scale random intercept and slope
scale_u_cov <- sqrt(1) * sqrt(0.09) * scale_rho
scale_u_sigma <- matrix(c(1, scale_u_cov,
                          scale_u_cov, 0.09), nrow = 2)

# Generate location random effects
location_raneff <- mvrnorm(Nnum, mu = location_beta, Sigma = location_u_sigma)
location_raneff <- as.data.frame(location_raneff)
names(location_raneff) <- c("mu_i", "phi_i")
location_raneff$ID <- 1:Nnum

# Generate scale random effects
scale_raneff <- mvrnorm(Nnum, mu = scale_beta, Sigma = scale_u_sigma)
scale_raneff <- as.data.frame(scale_raneff)
names(scale_raneff) <- c("alpha_i", "beta_i")
scale_raneff$ID <- 1:Nnum

# Create longitudinal dataset
mydata <- data.frame(matrix(NA, nrow = Nnum * Tnum, ncol = 2))
names(mydata) <- c("ID", "t")
mydata$ID <- rep(1:Nnum, each = Tnum)
mydata$t <- rep(1:Tnum, Nnum)

# Generate predictor and error terms
mydata$x_it <- rnorm(n = Nnum * Tnum, mean = 0, sd = 1)

# Merge random effects with main data
mydata <- left_join(mydata, location_raneff[, c("ID", "mu_i", "phi_i")], by = "ID")
mydata <- left_join(mydata, scale_raneff[, c("ID", "alpha_i", "beta_i")], by = "ID")

# Generate outcome with heteroscedastic errors
for(i in 1:nrow(mydata)) {
  # Compute conditional variance: exp(alpha_i + beta_i * x_it)
  conditional_var <- exp(mydata$alpha_i[i] + mydata$beta_i[i] * mydata$x_it[i])
  # Generate heteroscedastic error
  mydata$e_it[i] <- rnorm(1, mean = 0, sd = sqrt(conditional_var))
}

# Generate outcome: location component + heteroscedastic error
mydata$y_it <- mydata$mu_i + mydata$phi_i * mydata$x_it + mydata$e_it

# Summary statistics
summary(location_raneff)
summary(scale_raneff)
summary(mydata)

## 2 Parameter Estimation in JAGS ##############################################

# Load required packages
library(rjags)

jags_text <- "
model {
  # Likelihood ---------------------------------
  for (j in 1:Nobs) {
    y[j] ~ dnorm(y_mu[j], y_pre[j])
    y_mu[j] <- mu_location[subject[j]] + phi_location[subject[j]] * x[j]
    y_pre[j] <- exp(-(alpha_scale[subject[j]] + beta_scale[subject[j]] * x[j]))
  }

  # Location random effects --------------------
  for (i in 1:N) {
    mu_location[i] <- location_raneff[i,1]
    phi_location[i] <- location_raneff[i,2]
    location_raneff[i,1:2] ~ dmnorm(location_beta[1:2], location_pre[1:2,1:2])
  }

  # Scale random effects -----------------------
  for (i in 1:N) {
    alpha_scale[i] <- scale_raneff[i,1]
    beta_scale[i] <- scale_raneff[i,2]
    scale_raneff[i,1:2] ~ dmnorm(scale_beta[1:2], scale_pre[1:2,1:2])
  }

  # Priors for location parameters ------------
  location_beta[1] ~ dnorm(0, 0.01)    # Location intercept fixed effect
  location_beta[2] ~ dnorm(0, 0.01)    # Location slope fixed effect

  location_tau_beta_0 ~ dgamma(0.001, 0.001)
  location_sd_beta_0 <- sqrt(1/location_tau_beta_0)
  location_tau_beta_1 ~ dgamma(0.001, 0.001)
  location_sd_beta_1 <- sqrt(1/location_tau_beta_1)
  location_rho ~ dunif(-1, 1)

  location_sigma[1,1] <- location_sd_beta_0 * location_sd_beta_0
  location_sigma[2,2] <- location_sd_beta_1 * location_sd_beta_1
  location_sigma[1,2] <- location_sd_beta_0 * location_rho * location_sd_beta_1
  location_sigma[2,1] <- location_sigma[1,2]
  location_pre[1:2,1:2] <- inverse(location_sigma[1:2,1:2])

  # Priors for scale parameters --------------
  scale_beta[1] ~ dnorm(0, 0.01)      # Scale intercept fixed effect
  scale_beta[2] ~ dnorm(0, 0.01)      # Scale slope fixed effect

  scale_tau_beta_0 ~ dgamma(0.001, 0.001)
  scale_sd_beta_0 <- sqrt(1/scale_tau_beta_0)
  scale_tau_beta_1 ~ dgamma(0.001, 0.001)
  scale_sd_beta_1 <- sqrt(1/scale_tau_beta_1)
  scale_rho ~ dunif(-1, 1)

  scale_sigma[1,1] <- scale_sd_beta_0 * scale_sd_beta_0
  scale_sigma[2,2] <- scale_sd_beta_1 * scale_sd_beta_1
  scale_sigma[1,2] <- scale_sd_beta_0 * scale_rho * scale_sd_beta_1
  scale_sigma[2,1] <- scale_sigma[1,2]
  scale_pre[1:2,1:2] <- inverse(scale_sigma[1:2,1:2])
}
"

data_list <- list(
  Nobs = nrow(mydata),
  subject = mydata$ID,
  N = length(unique(mydata$ID)),
  y = mydata$y_it,
  x = mydata$x_it
)

# Initialize the JAGS model
jags_file <- textConnection(jags_text)
jags_model <- jags.model(jags_file,
                         data = data_list,
                         n.chain = 4)
# Run the MCMC sampler with a burn-in period of 5000 iterations
update(jags_model, 5000)
# Draw posterior samples
variables <- c("mu_location", "phi_location", "alpha_scale", "beta_scale",
               "location_beta", "scale_beta",
               "location_tau_beta_0", "location_tau_beta_1", "location_rho",
               "scale_tau_beta_0", "scale_tau_beta_1", "scale_rho")
post <- coda.samples(jags_model,
                     variable.names = variables,
                     n.iter = 25000, thin = 10)

# Combine posterior samples from all chains into a single matrix
samps <- do.call(rbind, post)
dim(samps)
##| 10000   810

# Extract posterior samples for model parameters only (excluding latent variables)
# 10 model parameters
pars_vector <- c("location_beta[1]", "location_beta[2]",
                 "scale_beta[1]", "scale_beta[2]",
                 "location_tau_beta_0", "location_tau_beta_1", "location_rho",
                 "scale_tau_beta_0", "scale_tau_beta_1", "scale_rho")
samps2 <- as.matrix(samps[, pars_vector])
dim(samps2)
##| 10000   10

## 3 Model Evaluation with bleval ##############################################

# Install bleval if not already installed
# devtools::install_github("luoxh3/bleval")
library(bleval)

## 3.1 Compute Information Criteria --------------------------------------------

## Step 1: Specify the log_joint_i function

log_joint_i <- function(samples_s, data, i, Ngrid, nodes) {
  # ==================================================
  # STEP 1: Extract data for the current unit i
  # ==================================================
  Nobs <- data$Nobs    # Total number of observations
  Nnum <- data$N       # Number of level-2 units
  Tnum <- Nobs / Nnum  # Number of observations per unit

  # Extract observations for unit i
  y_i <- data$y[((i-1)*Tnum+1):(i*Tnum)]  # Outcome values for unit i
  x_i <- data$x[((i-1)*Tnum+1):(i*Tnum)]  # Predictor values for unit i

  # ==================================================
  # STEP 2: Prepare data structures for quadrature
  # ==================================================
  # The 'nodes' matrix contains quadrature points for 4 latent variables
  # Each row represents one combination of latent variable values
  # For 4 latent variables with Ngrid = 5: nodes has 625 rows, 4 columns

  # Expand x and y to match the quadrature grid dimensions
  total_nodes <- nrow(nodes)
  x_i_extended_mat <- matrix(rep(x_i, times = total_nodes),
                             nrow = total_nodes, byrow = TRUE)
  y_i_extended_mat <- matrix(rep(y_i, times = total_nodes),
                             nrow = total_nodes, byrow = TRUE)

  # ==================================================
  # STEP 3: Compute log conditional likelihood
  # ==================================================
  # This corresponds to: log p(y_j | η_j, θ)
  # In location-scale model:
  # y_it ~ N(μ_location_i + φ_location_i * x_it, exp(α_scale_i + β_scale_i * x_it))

  # Extract latent variables from nodes matrix
  mu_location <- nodes[, 1]    # Location random intercept
  phi_location <- nodes[, 2]   # Location random slope
  alpha_scale <- nodes[, 3]    # Scale random intercept
  beta_scale <- nodes[, 4]     # Scale random slope

  # Compute predicted means and variances for each quadrature point
  predicted_mean <- mu_location + phi_location * x_i_extended_mat
  predicted_var <- exp(alpha_scale + beta_scale * x_i_extended_mat)

  # Compute log-density for each observation at each quadrature point
  log_con_t_i <- dnorm(y_i_extended_mat, mean = predicted_mean,
                       sd = sqrt(predicted_var), log = TRUE)

  # Sum over observations within the unit (conditional independence)
  log_con_i <- rowSums(log_con_t_i)

  # ==================================================
  # STEP 4: Compute log density for latent variables
  # ==================================================
  # This corresponds to: log p(η_j | θ)
  # We have two independent bivariate normal distributions:
  # location_raneff ~ MVN(location_beta, location_sigma)
  # scale_raneff ~ MVN(scale_beta, scale_sigma)

  # Location random effects distribution
  location_sd_0 <- sqrt(1 / samples_s[["location_tau_beta_0"]])
  location_sd_1 <- sqrt(1 / samples_s[["location_tau_beta_1"]])
  location_mean <- c(samples_s[["location_beta[1]"]], samples_s[["location_beta[2]"]])

  location_sigma <- matrix(c(
    location_sd_0^2,
    location_sd_0 * samples_s[["location_rho"]] * location_sd_1,
    location_sd_0 * samples_s[["location_rho"]] * location_sd_1,
    location_sd_1^2
  ), nrow = 2)

  # Scale random effects distribution
  scale_sd_0 <- sqrt(1 / samples_s[["scale_tau_beta_0"]])
  scale_sd_1 <- sqrt(1 / samples_s[["scale_tau_beta_1"]])
  scale_mean <- c(samples_s[["scale_beta[1]"]], samples_s[["scale_beta[2]"]])

  scale_sigma <- matrix(c(
    scale_sd_0^2,
    scale_sd_0 * samples_s[["scale_rho"]] * scale_sd_1,
    scale_sd_0 * samples_s[["scale_rho"]] * scale_sd_1,
    scale_sd_1^2
  ), nrow = 2)

  # Extract location and scale components from nodes
  location_nodes <- nodes[, 1:2]
  scale_nodes <- nodes[, 3:4]

  # Compute multivariate normal log-densities
  log_location_prior <- mvtnorm::dmvnorm(location_nodes, location_mean, location_sigma, log = TRUE)
  log_scale_prior <- mvtnorm::dmvnorm(scale_nodes, scale_mean, scale_sigma, log = TRUE)

  # Total latent variable prior (independent components)
  log_raneff_i <- log_location_prior + log_scale_prior

  # ==================================================
  # STEP 5: Return the log joint density for each unit
  # ==================================================
  # This equals: log p(y_j | η_j, θ) + log p(η_j | θ)
  log_con_i + log_raneff_i
}

## Step 2: Compute the posterior means and covariance matrices of latent variables
# We have 4 latent variables per unit: mu_location, phi_location, alpha_scale, beta_scale
raneff_mu_list <- vector("list", Nnum)
raneff_cov_list <- vector("list", Nnum)

for (i in 1:Nnum) {
  # Extract posterior means for all 4 random effects
  raneff_mu_list[[i]] <- c(
    mean(samps[, paste0("mu_location[", i, "]")]),
    mean(samps[, paste0("phi_location[", i, "]")]),
    mean(samps[, paste0("alpha_scale[", i, "]")]),
    mean(samps[, paste0("beta_scale[", i, "]")])
  )

  # Extract posterior covariance matrix for all 4 random effects
  raneff_cov_list[[i]] <- cov(samps[, c(
    paste0("mu_location[", i, "]"),
    paste0("phi_location[", i, "]"),
    paste0("alpha_scale[", i, "]"),
    paste0("beta_scale[", i, "]")
  )])
}

## Step 3: Compute the log marginal likelihood
# Note: Using smaller Ngrid (=5) due to 4-dimensional integration
log_marglik_result <- bleval::log_marglik(
  samples = samps2,
  data = data_list,
  Ngrid = 5,
  lv_mu = raneff_mu_list,
  lv_cov = raneff_cov_list,
  log_joint_i = log_joint_i
)

## Step 4: Compute information criteria
IC_results <- bleval::calc_IC(log_marglik_result, 1)
print(IC_results)
##|    p_dic      elpd_dic           dic        p_waic     elpd_waic          waic       p_looic    elpd_looic         looic
##| 9.992118 -15103.784887  30207.569774     10.518208 -15104.170101  30208.340201     10.530743 -15104.182635  30208.365271

## 3.2 Compute Fully Marginal Likelihood ---------------------------------------

## Step 1: Specify the log_prior function
log_prior <- function(samples_s) {
  dnorm(samples_s[["location_beta[1]"]], mean = 0, sd = sqrt(1/0.01), log = TRUE) +
  dnorm(samples_s[["location_beta[2]"]], mean = 0, sd = sqrt(1/0.01), log = TRUE) +
  dnorm(samples_s[["scale_beta[1]"]], mean = 0, sd = sqrt(1/0.01), log = TRUE) +
  dnorm(samples_s[["scale_beta[2]"]], mean = 0, sd = sqrt(1/0.01), log = TRUE) +
  dgamma(samples_s[["location_tau_beta_0"]], shape = 0.001, rate = 0.001, log = TRUE) +
  dgamma(samples_s[["location_tau_beta_1"]], shape = 0.001, rate = 0.001, log = TRUE) +
  dgamma(samples_s[["scale_tau_beta_0"]], shape = 0.001, rate = 0.001, log = TRUE) +
  dgamma(samples_s[["scale_tau_beta_1"]], shape = 0.001, rate = 0.001, log = TRUE) +
  dunif(samples_s[["location_rho"]], min = -1, max = 1, log = TRUE) +
  dunif(samples_s[["scale_rho"]], min = -1, max = 1, log = TRUE)
}

## Step 2: Define parameter bounds
lb <- c(rep(-Inf, 4), rep(0, 2), -1, rep(0, 2), -1)
ub <- c(rep(Inf, 4), rep(Inf, 2), 1, rep(Inf, 2), 1)
names(lb) <- pars_vector
names(ub) <- pars_vector

## Step 3: Compute the log fully marginal likelihood
log_fmarglik_result <- bleval::log_fmarglik(
  samples = samps2,
  data = data_list,
  Ngrid = 5,
  lv_mu = raneff_mu_list,
  lv_cov = raneff_cov_list,
  log_joint_i = log_joint_i,
  log_prior = log_prior,
  lb = lb,
  ub = ub
)

print(log_fmarglik_result)
##| Bridge sampling estimate of the log marginal likelihood: -15152.77
##| Estimate obtained in 4 iteration(s) via method "normal".
