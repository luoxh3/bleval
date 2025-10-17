## This R script demonstrates Bayesian evaluation of a latent mediation model.
##
## Key characteristics:
## - Model category: Structural Equation Model (SEM)
## - Number of latent variables: 3

## 1 Data Generation ###########################################################

# Load required packages
library(lavaan)

set.seed(123)

N <- 240  # Number of respondents

# Specify the data-generating latent mediation model:
# For X (x1-x3):
# - item factor loadings: c(1.00, 0.80, 1.20)
# - item error variances: c(0.36, 0.36, 0.36)
# - factor variance: 1.00
# - mu_X = 0
# For M (x4-x6):
# - item factor loadings: c(1.00, 1.20, 0.80)
# - item error variances: c(0.36, 0.36, 0.49)
# - factor residual variance: 1.00
# - mu_M = 0.8 * X
# For Y (x7-x10):
# - item factor loadings: c(1.00, 0.80, 1.00, 1.20)
# - item error variances: c(0.49, 0.49, 0.49, 0.49)
# - factor residual variance: 1.00
# - mu_Y = 0.2 * X + 0.5 * M
population_model <- ' X =~ x1 + 0.8*x2 + 1.2*x3
                      M =~ m1 + 1.2*m2 + 0.8*m3
                      Y =~ y1 + 0.8*y2 + y3 + 1.2*y4
                      X ~~ 1*X
                      M ~~ 1*M
                      Y ~~ 1*Y
                      x1 ~~ 0.36*x1
                      x2 ~~ 0.36*x2
                      x3 ~~ 0.36*x3
                      m1 ~~ 0.36*m1
                      m2 ~~ 0.36*m2
                      m3 ~~ 0.49*m3
                      y1 ~~ 0.49*y1
                      y2 ~~ 0.49*y2
                      y3 ~~ 0.49*y3
                      y4 ~~ 0.49*y4
                      Y ~ 0.2 * X + 0.5 * M
                      M ~ 0.8 * X
                    '
# Generate CFA dataset
data <- simulateData(
  model = population_model, 
  sample.nobs = N,
  meanstructure = FALSE
)

# Summary statistics
summary(data)

## 2 Parameter Estimation in JAGS ##############################################

# Load required packages
library(rjags)

jags_text <- "
model {
  # likelihood ---------------------------------
  # Measurement Models
  for(i in 1:N) {
    eta_X[i] ~ dnorm(0, tau_eta_X)
    eta_M[i] ~ dnorm(mu_eta_M[i], tau_eta_M)
    eta_Y[i] ~ dnorm(mu_eta_Y[i], tau_eta_Y)
  
    for(j in 1:NumX) {
      X[i,j] ~ dnorm(mu_X[i,j], tau_X[j])
      mu_X[i,j] <- lambda_X[j] * eta_X[i]
    }
    for(j in 1:NumM) {
      M[i,j] ~ dnorm(mu_M[i,j], tau_M[j])
      mu_M[i,j] <- lambda_M[j] * eta_M[i]
    }
    for(j in 1:NumY) {
      Y[i,j] ~ dnorm(mu_Y[i,j], tau_Y[j])
      mu_Y[i,j] <- lambda_Y[j] * eta_Y[i]
    }
  }
  # Structural Model
  for (i in 1:N) {
    mu_eta_M[i] <- a * eta_X[i]
    mu_eta_Y[i] <- c_prime * eta_X[i] + b * eta_M[i]
  }
  
  # Priors -------------------------------------
  a ~ dnorm(0, 0.01)
  b ~ dnorm(0, 0.01)
  c_prime ~ dnorm(0, 0.01)
  
  lambda_X[1] <- 1
  lambda_M[1] <- 1
  lambda_Y[1] <- 1
  for(j in 2:NumX) {
    lambda_X[j] ~ dnorm(0, 0.01)
  }
  for(j in 2:NumM) {
    lambda_M[j] ~ dnorm(0, 0.01)
  }
  for(j in 2:NumY) {
    lambda_Y[j] ~ dnorm(0, 0.01)
  }
  
  for(j in 1:NumX) {
    tau_X[j] ~ dgamma(0.001, 0.001)
  }
  for(j in 1:NumM) {
    tau_M[j] ~ dgamma(0.001, 0.001)
  }
  for(j in 1:NumY) {
    tau_Y[j] ~ dgamma(0.001, 0.001)
  }
  
  tau_eta_X ~ dgamma(0.001, 0.001)
  tau_eta_M ~ dgamma(0.001, 0.001)
  tau_eta_Y ~ dgamma(0.001, 0.001)
}"

X <- as.matrix(data[, c("x1", "x2", "x3")])
M <- as.matrix(data[, c("m1", "m2", "m3")])
Y <- as.matrix(data[, c("y1", "y2", "y3", "y4")])

NumX <- ncol(X)
NumM <- ncol(M)
NumY <- ncol(Y)

data_list <- list(
  X = X,
  M = M,
  Y = Y,
  NumX = NumX,
  NumM = NumM,
  NumY = NumY,
  N = N
)

# Initialize the JAGS model
jags_file <- textConnection(jags_text)
jags_model <- jags.model(jags_file,
                         data = data_list,
                         n.chain = 4)
# Run the MCMC sampler with a burn-in period of 5000 iterations
update(jags_model, 5000)
# Draw posterior samples
post <- coda.samples(jags_model,
                     variable.names = c("a", "b", "c_prime", 
                                        "lambda_X", "lambda_M", "lambda_Y",
                                        "tau_X", "tau_M", "tau_Y",
                                        "tau_eta_X", "tau_eta_M", "tau_eta_Y",
                                        "eta_X", "eta_M", "eta_Y"),
                     n.iter = 25000, thin = 10)

# Combine posterior samples from all chains into a single matrix
samps <- do.call(rbind, post)

# Extract posterior samples for model parameters only (excluding latent variables)
# 23 model parameters
pars_vector <- c("a", "b", "c_prime",
                 "lambda_X[2]","lambda_X[3]","lambda_M[2]","lambda_M[3]",
                 "lambda_Y[2]","lambda_Y[3]","lambda_Y[4]",
                 "tau_X[1]", "tau_X[2]","tau_X[3]",
                 "tau_M[1]","tau_M[2]","tau_M[3]",
                 "tau_Y[1]","tau_Y[2]","tau_Y[3]","tau_Y[4]",
                 "tau_eta_X", "tau_eta_M", "tau_eta_Y")
samps2 <- as.matrix(samps[, pars_vector])
dim(samps2)
##| 10000   23

## 3 Model Evaluation with bleval ##############################################

# Install bleval if not already installed
# devtools::install_github("luoxh3/bleval")
library(bleval)

## 3.1 Compute Information Criteria --------------------------------------------

## Step 1: Specify the log_joint_i function

log_joint_i <- function(samples_s, data, i, Ngrid, nodes) {
  # ==================================================
  # STEP 1: Extract data and parameter draws needed
  # ==================================================
  X_i <- data$X[i,]
  M_i <- data$M[i,]
  Y_i <- data$Y[i,]
  
  Ndim <- 3

  # Extract posterior draws for parameters
  a <- samples_s[["a"]]
  b <- samples_s[["b"]]
  c_prime <- samples_s[["c_prime"]]
  
  lambda_X <- as.vector(samples_s[grep("^lambda_X\\[\\d+\\]$", names(samples_s))])
  lambda_X <- c(1, lambda_X)
  lambda_M <- as.vector(samples_s[grep("^lambda_M\\[\\d+\\]$", names(samples_s))])
  lambda_M <- c(1, lambda_M)
  lambda_Y <- as.vector(samples_s[grep("^lambda_Y\\[\\d+\\]$", names(samples_s))])
  lambda_Y <- c(1, lambda_Y)
  
  sd_X_mat <- matrix(rep( sqrt(1/samples_s[grep("^tau_X\\[\\d+\\]$", names(samples_s))]), 
                          times = Ngrid^Ndim), nrow = Ngrid^Ndim, byrow = TRUE)
  sd_M_mat <- matrix(rep( sqrt(1/samples_s[grep("^tau_M\\[\\d+\\]$", names(samples_s))]), 
                          times = Ngrid^Ndim), nrow = Ngrid^Ndim, byrow = TRUE)
  sd_Y_mat <- matrix(rep( sqrt(1/samples_s[grep("^tau_Y\\[\\d+\\]$", names(samples_s))]), 
                          times = Ngrid^Ndim), nrow = Ngrid^Ndim, byrow = TRUE)
  
  sd_eta_X <- sqrt(1/samples_s[["tau_eta_X"]])
  sd_eta_M <- sqrt(1/samples_s[["tau_eta_M"]])
  sd_eta_Y <- sqrt(1/samples_s[["tau_eta_Y"]])
  
  # ==================================================
  # STEP 2: Prepare data structures for quadrature
  # ==================================================
  # Dimensions of nodes: Ngrid^Ndim rows, Ndim columns (Ndim = number of latent variables)
  # Ndim = 3
  
  # Expand the data to match the quadrature grid dimensions
  X_i_extended_mat <- matrix(rep(X_i, times = Ngrid^Ndim), 
                             nrow = Ngrid^Ndim, ncol = data$NumX, byrow = TRUE)
  M_i_extended_mat <- matrix(rep(M_i, times = Ngrid^Ndim), 
                             nrow = Ngrid^Ndim, ncol = data$NumM, byrow = TRUE)
  Y_i_extended_mat <- matrix(rep(Y_i, times = Ngrid^Ndim), 
                             nrow = Ngrid^Ndim, ncol = data$NumY, byrow = TRUE)
  
  # ==================================================
  # STEP 3: Compute log conditional likelihood
  # ==================================================
  # This corresponds to: log p(y_j | η_j, θ)
  
  # Compute log-density for each observation at each quadrature point
  log_con_X_ij <- dnorm(X_i_extended_mat, mean = nodes[,1] %*% t(lambda_X), 
                         sd = sd_X_mat, log = TRUE)
  log_con_M_ij <- dnorm(M_i_extended_mat, mean = nodes[,2] %*% t(lambda_M), 
                         sd = sd_M_mat, log = TRUE)
  log_con_Y_ij <- dnorm(Y_i_extended_mat, mean = nodes[,3] %*% t(lambda_Y), 
                         sd = sd_Y_mat, log = TRUE)
  
  # Sum over observations within the unit (conditional independence)
  log_con_i <- rowSums(log_con_X_ij) + rowSums(log_con_M_ij) + rowSums(log_con_Y_ij)
  
  # ==================================================
  # STEP 4: Compute log density for latent variables
  # ==================================================
  # This corresponds to: log p(η_j | θ)
  predicted_M <- samples_s[["a"]] * nodes[,1] 
  predicted_Y <- samples_s[["c_prime"]] * nodes[,1] + samples_s[["b"]] * nodes[,2]
  log_eta_X_i <- dnorm(nodes[,1], mean = 0, sd = sd_eta_X, log = TRUE)
  log_eta_M_i <- dnorm(nodes[,2], mean = predicted_M, sd = sd_eta_M, log = TRUE)
  log_eta_Y_i <- dnorm(nodes[,3], mean = predicted_Y, sd = sd_eta_Y, log = TRUE)
  
  # ==================================================
  # STEP 5: Return the log joint density for each unit
  # ==================================================
  # This equals: log p(y_j | η_j, θ) + log p(η_j | θ)
  log_con_i + (log_eta_X_i + log_eta_M_i + log_eta_Y_i)
}

## Step 2: Compute the posterior means and covariance matrices of latent variables
# Create lists to store posterior draws of latent variable (eta)
lv_list <- list()
for (i in 1:N){
  lv_list[[i]] <- cbind(samps[, paste0("eta_X[", i, "]")],
                        samps[, paste0("eta_M[", i, "]")],
                        samps[, paste0("eta_Y[", i, "]")])
}
# Compute the posterior means and covariance matrices
lv_mu_list <- lapply(lv_list, colMeans)
lv_cov_list <- lapply(lv_list, cov)

## Step 3: Compute the log marginal likelihood
# Note: Using smaller Ngrid (=5) due to 3-dimensional integration
log_marglik_result <- bleval::log_marglik(
  samples = samps2,
  data = data_list, 
  Ngrid = 5,
  lv_mu = lv_mu_list, 
  lv_cov = lv_cov_list,
  log_joint_i = log_joint_i
)

## Step 4: Compute information criteria
IC_results <- bleval::calc_IC(log_marglik_result, 1)
print(IC_results)
##|    p_dic    elpd_dic         dic      p_waic   elpd_waic        waic     p_looic  elpd_looic       looic 
##| 23.20614 -3162.88165  6325.76329    22.86467 -3162.93871  6325.87742    22.89046 -3162.96450  6325.92899

## 3.2 Compute Fully Marginal Likelihood ---------------------------------------

## Step 1: Specify the log_prior function
log_prior <- function(samples_s) {
  a <- samples_s[["a"]]
  b <- samples_s[["b"]]
  c_prime <- samples_s[["c_prime"]]
  
  lambda_X <- as.vector(samples_s[grep("^lambda_X\\[\\d+\\]$", names(samples_s))])
  lambda_M <- as.vector(samples_s[grep("^lambda_M\\[\\d+\\]$", names(samples_s))])
  lambda_Y <- as.vector(samples_s[grep("^lambda_Y\\[\\d+\\]$", names(samples_s))])
  
  tau_X <- as.vector(samples_s[grep("^tau_X\\[\\d+\\]$", names(samples_s))])
  tau_M <- as.vector(samples_s[grep("^tau_M\\[\\d+\\]$", names(samples_s))])
  tau_Y <- as.vector(samples_s[grep("^tau_Y\\[\\d+\\]$", names(samples_s))])

  tau_eta_X <- samples_s[["tau_eta_X"]]
  tau_eta_M <- samples_s[["tau_eta_M"]]
  tau_eta_Y <- samples_s[["tau_eta_Y"]]
  
  dnorm(a, mean = 0, sd = sqrt(1/0.01), log = TRUE) +
    dnorm(b, mean = 0, sd = sqrt(1/0.01), log = TRUE) +
    dnorm(c_prime, mean = 0, sd = sqrt(1/0.01), log = TRUE) +
    sum(dnorm(lambda_X, mean = 0, sd = sqrt(1/0.01), log = TRUE)) +
    sum(dnorm(lambda_M, mean = 0, sd = sqrt(1/0.01), log = TRUE)) +
    sum(dnorm(lambda_Y, mean = 0, sd = sqrt(1/0.01), log = TRUE)) +
    sum(dgamma(tau_X, shape = 0.001, rate = 0.001, log = TRUE)) +
    sum(dgamma(tau_M, shape = 0.001, rate = 0.001, log = TRUE)) +
    sum(dgamma(tau_Y, shape = 0.001, rate = 0.001, log = TRUE)) +
    dgamma(tau_eta_X, shape = 0.001, rate = 0.001, log = TRUE) +
    dgamma(tau_eta_M, shape = 0.001, rate = 0.001, log = TRUE) +
    dgamma(tau_eta_Y, shape = 0.001, rate = 0.001, log = TRUE)
}

## Step 2: Define parameter bounds
lb <- c(rep(-Inf, 10), rep(0, 13))
ub <- c(rep(Inf, 10), rep(Inf, 13))
names(lb) <- pars_vector
names(ub) <- pars_vector

## Step 3: Compute the log fully marginal likelihood
log_fmarglik_result <- bleval::log_fmarglik(
  samples = samps2, 
  data = data_list, 
  Ngrid = 5,
  lv_mu = lv_mu_list, 
  lv_cov = lv_cov_list,
  log_joint_i = log_joint_i, 
  log_prior = log_prior,
  lb = lb,
  ub = ub
)

print(log_fmarglik_result)
##| Bridge sampling estimate of the log marginal likelihood: -3295.84
##| Estimate obtained in 6 iteration(s) via method "normal".
