## This R script demonstrates Bayesian evaluation of a 2-Parameter Logistic (2-PL) model.
##
## Key characteristics:
## - Model category: Item Response Theory (IRT) model
## - Number of latent variables: 1

## 1 Data Generation ###########################################################

# Load required packages
library(mirt)

set.seed(123)

N <- 300  # Number of respondents
J <- 6   # Number of items

# Specify the data-generating 2-PL model
a <- c(1.0, 1.5, 2.0, 1.0, 1.5, 2.0)  # item discrimination parameters
b <- c(-0.5, -0.5, 0, 0, 0.5, 0.5) # item difficulty parameters

# Generate 2-PL dataset
data <- simdata(a = a, d = -a*b, N = N, itemtype = '2PL')

# Summary statistics
summary(data)

## 2 Parameter Estimation in JAGS ##############################################

# Load required packages
library(rjags)

jags_text <- "
model {
  # likelihood ---------------------------------
  for (i in 1:N) {
    theta[i] ~ dnorm(0, 1)
    
    for (j in 1:J) {
      logit(p[i,j]) <- a[j] * (theta[i] - b[j])
      y[i,j] ~ dbern(p[i,j])
    }
  }
  
  # Priors -------------------------------------
  for (j in 1:J) {
    a[j] ~ dgamma(0.001, 0.001)
    b[j] ~ dnorm(0, 0.01)
  }
}"

data_list <- list(
  N = N,
  J = J,
  y = data
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
                     variable.names = c("a", "b", "theta"),
                     n.iter = 25000, thin = 10)

# Combine posterior samples from all chains into a single matrix
samps <- do.call(rbind, post)

# Extract posterior samples for model parameters only (excluding latent variables)
# 12 model parameters
pars_vector <- c("a[1]", "a[2]", "a[3]", "a[4]", "a[5]", "a[6]",
                 "b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]")
samps2 <- as.matrix(samps[, pars_vector])
dim(samps2)
##| 10000   12

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
  a <- as.vector(samples_s[grep("^a\\[\\d+\\]$", names(samples_s))])
  b <- as.vector(samples_s[grep("^b\\[\\d+\\]$", names(samples_s))])
  
  # ==================================================
  # STEP 2: Prepare data structures for quadrature
  # ==================================================
  # Dimensions of nodes: Ngrid^Ndim rows, Ndim columns (Ndim = number of latent variables)
  # Since the model contains only 1 latent variable, the 'nodes' is a column vector of quadrature points
  
  # Expand data and parameters to match the quadrature grid dimensions
  y_i_extended_mat <- matrix(rep(y_i, times = Ngrid^1),
                             nrow = Ngrid^1, ncol = J, byrow = TRUE)
  
  a_extended_mat <- matrix(rep(a, times = Ngrid^1),
                           nrow = Ngrid^1, ncol = J, byrow = TRUE)
  b_extended_mat <- matrix(rep(b, times = Ngrid^1),
                           nrow = Ngrid^1, ncol = J, byrow = TRUE)
  theta_i_extended_mat <- matrix(rep(nodes, each = J),
                                 nrow = Ngrid^1, ncol = J, byrow = TRUE)
  
  # ==================================================
  # STEP 3: Compute log conditional likelihood
  # ==================================================
  # This corresponds to: log p(y_j | η_j, θ)
  # In 2-PL model, for person i and item j, y_ij ~ Bernoulli(a_j*(theta_i - b_j))
  
  # Compute log-density for each observation at each quadrature point
  p <- plogis(a_extended_mat * (theta_i_extended_mat - b_extended_mat))
  log_con_i_j <- dbinom(y_i_extended_mat, 1, p, log = TRUE)
  
  # Sum over observations within the unit (conditional independence)
  log_con_i <- rowSums(log_con_i_j)
  
  # ==================================================
  # STEP 4: Compute log density for latent variables
  # ==================================================
  # This corresponds to: log p(η_j | θ)
  # In 2-PL, we have a normal distribution for person i's latent ability: theta_i ~ N(0, 1)
  
  # Compute log normal density
  log_lv_i <- dnorm(nodes, mean = 0, sd = 1, log = TRUE)
  
  # ==================================================
  # STEP 5: Return the log joint density for each unit
  # ==================================================
  # This equals: log p(y_j | η_j, θ) + log p(η_j | θ)
  log_con_i + log_lv_i
}

## Step 2: Compute the posterior means and covariance matrices of latent variables
# Create lists to store posterior draws of latent variable (eta)
lv_list <- list()
for (i in 1:N){
  lv_list[[i]] <- samps[, paste0("theta[", i, "]")]
}
# For the only latent variable (theta), compute the posterior mean and variance
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
##| 12.35940 -1091.90435  2183.80871    12.21365 -1091.90767  2183.81533    12.22524 -1091.91925  2183.83851  

## 3.2 Compute Fully Marginal Likelihood ---------------------------------------

## Step 1: Specify the log_prior function
log_prior <- function(samples_s) {
  a <- as.vector(samples_s[grep("^a\\[\\d+\\]$", names(samples_s))])
  b <- as.vector(samples_s[grep("^b\\[\\d+\\]$", names(samples_s))])
  
  sum(dnorm(b, mean = 0, sd = sqrt(1/0.01), log = TRUE)) +
  sum(dgamma(a, shape = 0.001, rate = 0.001, log = TRUE))
}

## Step 2: Define parameter bounds
lb <- c(rep(0, J), rep(-Inf, J))
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
##| Bridge sampling estimate of the log marginal likelihood: -1152.107
##| Estimate obtained in 5 iteration(s) via method "normal".
