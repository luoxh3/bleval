## This R script demonstrates Bayesian evaluation of a confirmatory factor analysis 
## (CFA) model.
##
## Key characteristics:
## - Model category: Structural Equation Model (SEM)
## - Number of latent variables: 1

## 1 Data Generation ###########################################################

# Load required packages
library(lavaan)

set.seed(123)

N <- 100  # Number of respondents
J <- 4   # Number of items

# Specify the data-generating CFA model:
# - item factor loadings: c(1.00, 1.20, 1.00, 0.80)
# - item error variances: c(0.36, 0.36, 0.36, 0.36)
# - factor variance: 1.00
# - no mean structure
population_model <- ' f1 =~ x1 + 1.2*x2 + x3 + 0.8*x4
                      f1 ~~ 1*f1
                      x1 ~~ 0.36*x1
                      x2 ~~ 0.36*x2
                      x3 ~~ 0.36*x3
                      x4 ~~ 0.36*x4
                    '
# Generate CFA dataset
cfa_data <- simulateData(
  model = population_model, 
  sample.nobs = N,
  meanstructure = FALSE
)

# Summary statistics
summary(cfa_data)

## 2 Parameter Estimation in JAGS ##############################################

# Load required packages
library(rjags)

jags_text <- "
model {
  # likelihood ---------------------------------
   for(i in 1:N) {
    eta[i] ~ dnorm(0, tau_psi)
    
    for(j in 1:J) {
      y[i,j] ~ dnorm(mu[i,j], tau_theta[j])
      mu[i,j] <- lambda[j] * eta[i]
    }
   }
  
  # Priors -------------------------------------
  lambda[1] <- 1
  for(j in 2:J) {
    lambda[j] ~ dnorm(0, 0.01)
  }
  tau_psi ~ dgamma(0.001, 0.001)
  for(j in 1:J) {
    tau_theta[j] ~ dgamma(0.001, 0.001)
  }
}"

data_list <- list(
  N = N,
  J = J,
  y = cfa_data
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
                     variable.names = c("eta", "tau_psi", "lambda", "tau_theta"),
                     n.iter = 25000, thin = 10)

# Combine posterior samples from all chains into a single matrix
samps <- do.call(rbind, post)

# Extract posterior samples for model parameters only (excluding latent variables)
# 8 model parameters
pars_vector <- c("tau_psi", "lambda[2]", "lambda[3]", "lambda[4]",
                 "tau_theta[1]", "tau_theta[2]", "tau_theta[3]", "tau_theta[4]")
samps2 <- as.matrix(samps[, pars_vector])
dim(samps2)
##| 10000   8

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
  N <- data$N       # Number of respondents
  J <- data$J       # Number of items
  
  # Extract observations for unit i
  y <- data$y
  y_i <- as.matrix(y[i, ])

  # Extract posterior draws for parameters
  lambda <- as.vector(samples_s[grep("^lambda\\[\\d+\\]$", names(samples_s))])
  lambda <- c(1, lambda)
  tau_theta <- as.vector(samples_s[grep("^tau_theta\\[\\d+\\]$", names(samples_s))])
  tau_psi <- samples_s[["tau_psi"]]
  
  # ==================================================
  # STEP 2: Prepare data structures for quadrature
  # ==================================================
  # Dimensions of nodes: Ngrid^Ndim rows, Ndim columns (Ndim = number of latent variables)
  # Since the model contains only 1 latent variable, the 'nodes' is a column vector of quadrature points
  
  # Expand data to match the quadrature grid dimensions
  total_nodes <- length(nodes)
  y_i_extended_mat <- matrix(rep(y_i, times = total_nodes), 
                             nrow = total_nodes, byrow = TRUE) # Dims: total_nodes rows, J columns
  
  # ==================================================
  # STEP 3: Compute log conditional likelihood
  # ==================================================
  # This corresponds to: log p(y_j | η_j, θ)
  # In CFA model, for person i and item j, y_ij ~ N(lambda_j * eta_i, theta_j)
  
  # Compute log-density for each observation at each quadrature point
  log_con_i_j <- dnorm(y_i_extended_mat, mean = nodes%*%t(lambda), 
                       sd = sqrt(rep(1/tau_theta, each = total_nodes)), log = TRUE)
  
  # Sum over observations within the unit (conditional independence)
  log_con_i <- rowSums(log_con_i_j)
  
  # ==================================================
  # STEP 4: Compute log density for latent variables
  # ==================================================
  # This corresponds to: log p(η_j | θ)
  # We have a normal distribution: eta_i ~ N(0, psi)
  
  # Compute log normal density
  log_lv_i <- dnorm(nodes, mean = 0, sd = sqrt(1/tau_psi), log = TRUE)
  
  # ==================================================
  # STEP 5: Return the log joint density for each unit
  # ==================================================
  # This equals: log p(y_j | η_j, θ) + log p(η_j | θ)
  log_con_i + log_lv_i
}

## Step 2: Compute the posterior means and covariance matrices of latent variables
# Create lists to store posterior draws of latent variable (eta)
lv_list <- list()
for (i in 1:N)
  lv_list[[i]] <- samps[, paste0("eta[", i, "]")]

# For the only latent variable (eta), compute the posterior mean and variance
lv_mu_list <- lapply(lv_list, mean)
lv_cov_list <- lapply(lv_list, var)

## Step 3: Compute the log marginal likelihood
log_marglik_result <- bleval::log_marglik(
  samples = samps2,
  data = data_list, 
  Ngrid = 9,
  lv_mu = lv_mu_list, 
  lv_cov = lv_cov_list,
  log_joint_i = log_joint_i
)

## Step 4: Compute information criteria
IC_results <- bleval::calc_IC(log_marglik_result, 1)
print(IC_results)
##|    p_dic    elpd_dic         dic      p_waic   elpd_waic        waic     p_looic  elpd_looic       looic 
##| 8.077615 -482.544284  965.088568    8.174065 -482.770144  965.540288    8.208325 -482.804404  965.608809 

## 3.2 Compute Fully Marginal Likelihood ---------------------------------------

## Step 1: Specify the log_prior function
log_prior <- function(samples_s) {
  lambda <- as.vector(samples_s[grep("^lambda\\[\\d+\\]$", names(samples_s))])
  tau_theta <- as.vector(samples_s[grep("^tau_theta\\[\\d+\\]$", names(samples_s))])
  tau_psi <- samples_s[["tau_psi"]]
  
  sum(dnorm(lambda, mean = 0, sd = sqrt(1/0.01), log = TRUE)) +
  sum(dgamma(tau_theta, shape = 0.001, rate = 0.001, log = TRUE)) +
  dgamma(tau_psi, shape = 0.001, rate = 0.001, log = TRUE)
}

## Step 2: Define parameter bounds
lb <- c(rep(-Inf, J), rep(0, J))
ub <- c(rep(Inf, J), rep(Inf, J))
names(lb) <- pars_vector
names(ub) <- pars_vector

## Step 3: Compute the log fully marginal likelihood
log_fmarglik_result <- bleval::log_fmarglik(
  samples = samps2, 
  data = data_list, 
  Ngrid = 9,
  lv_mu = lv_mu_list, 
  lv_cov = lv_cov_list,
  log_joint_i = log_joint_i, 
  log_prior = log_prior,
  lb = lb,
  ub = ub
)

print(log_fmarglik_result)
##| Bridge sampling estimate of the log marginal likelihood: -526.9813
##| Estimate obtained in 5 iteration(s) via method "normal".
