## This R script demonstrates Bayesian evaluation of a Graded Response Model (GRM).
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
K <- 5   # Number of categories

# Specify the data-generating GRM:
a <- c(1.0, 1.0, 1.5, 1.5, 2.0, 2.0) # Item discrimination parameters
b_1 <- c(-1.5, -1.0, -0.5, -1.5, -1.0, -0.5)
diff_mat <- matrix(runif(J * (K-2), 0.25, 0.75), nrow = J, ncol = K-2)
b_mat <- matrix(0, nrow = J, ncol = K-1)
b_mat[,1] <- b_1
for(k in 2:(K-1)) {
  b_mat[,k] <- b_mat[,k-1] + diff_mat[,k-1]
} # Item difficulty parameters

d_mat <- -a * b_mat

# Generate GRM dataset
data <- simdata(a = a, d = d_mat, N = N, itemtype = 'graded')

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
      y[i,j] ~ dcat(prob[i,j,1:K])
      
      for (k in 1:(K-1)) {
        logit(cumprob[i,j,k]) <- a[j] * (theta[i] - b[j,k])
      }
      
      prob[i,j,1] <- 1 - cumprob[i,j,1]                    # P(Y = 1)
      for (k in 2:(K-1)) {
        prob[i,j,k] <- cumprob[i,j,k-1] - cumprob[i,j,k]   # P(Y = k)
      }
      prob[i,j,K] <- cumprob[i,j,K-1]                      # P(Y = K)
    }
  }
  
  # Priors -------------------------------------
  for (j in 1:J) {
    a[j] ~ dlnorm(0, 1)
    b[j,1] ~ dnorm(0, 0.01)
    for (k in 2:(K-1)) {
      diff[j,k-1] ~ dnorm(0, 0.01)T(0,)
      b[j,k] <- b[j,k-1] + diff[j,k-1]
    }
  }
}"

data_list <- list(
  N = N,
  J = J,
  K = K,
  y = data + 1
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
                     variable.names = c("a", "b", "diff", "theta"),
                     n.iter = 25000, thin = 10)

# Combine posterior samples from all chains into a single matrix
samps <- do.call(rbind, post)

# Extract posterior samples for model parameters only (excluding latent variables)
# 30 model parameters
pars_vector <- c(paste("a[",1:J,"]", sep=""), paste("b[",1:J, ",1]", sep=""))
for (j in 1:J)
  for (k in 1:(K-2))
    pars_vector <- c(pars_vector, paste("diff[",j,",",k,"]", sep=""))

samps2 <- as.matrix(samps[, pars_vector])
dim(samps2)
##| 10000   30

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
  K <- data$K       # Number of categories
  
  # Extract observations for unit i
  y <- data$y
  y_i <- as.matrix(y[i, ])

  # Extract posterior draws for parameters
  a <- as.vector(samples_s[grep("^a\\[\\d+\\]$", names(samples_s))])
  b_1 <- as.vector(samples_s[grep("^b\\[\\d+,\\d+\\]$", names(samples_s))])
  diff <- matrix(samples_s[grep("^diff\\[\\d+,\\d+\\]$", names(samples_s))],
                 nrow = J, ncol = K-2, byrow = TRUE)
  b <- matrix(NA, nrow = J, ncol = K-1)
  b[,1] <- b_1
  for(k in 2:(K-1)) {
    b[,k] <- b[,k-1] + diff[,k-1]
  }
  
  # ==================================================
  # STEP 2: Prepare data structures for quadrature
  # ==================================================
  # Dimensions of nodes: Ngrid^Ndim rows, Ndim columns (Ndim = number of latent variables)
  # Since the model contains only 1 latent variable, the 'nodes' is a column vector of quadrature points
  
  # Note: 
  # STEP 2 was originally designed to expand data and parameters to enable efficient matrix computation in STEP 3. 
  # However, since the "prob" argument in the dcat() function only supports up to 2D arrays, 
  # we compute an Nnode x K probability matrix separately for each item j. 
  # This approach eliminates the need for pre-expanding data and parameters into Nnode x J matrices 
  # (e.g. y_i_extended_mat), and we instead process items sequentially for clarity.
  
  # ==================================================
  # STEP 3: Compute log conditional likelihood
  # ==================================================
  # This corresponds to: log p(y_j | η_j, θ)
  
  # For person i and item j, y_ij ~ Categorical(prob)
  # In GRM, category probabilities are computed from differences between adjacent cumulative probabilities
  
  prob <- matrix(NA, nrow = Ngrid^1, ncol = K)
  log_con_i_j <- matrix(NA, nrow = Ngrid^1, ncol = J)
  
  for (j in 1:J){
    # Compute category probabilities using GRM formulation
    prob[,1] <- 1 - plogis(a[j] * (nodes - b[j, 1]))
    for (k in 2:(K-1))
      prob[,k] <- plogis(a[j] * (nodes - b[j, k-1])) -  plogis(a[j] * (nodes - b[j, k]))
    prob[,K] <- plogis(a[j] * (nodes - b[j, K-1]))
    
    # Compute log-density for each observation at each quadrature point
    log_con_i_j[,j] <- extraDistr::dcat(y_i[j], prob, log = TRUE)
  }
  
  # Sum over observations within the unit (conditional independence)
  log_con_i <- rowSums(log_con_i_j)
  
  # ==================================================
  # STEP 4: Compute log density for latent variables
  # ==================================================
  # This corresponds to: log p(η_j | θ)
  
  # In GRM, we have a normal distribution for person i's latent trait: theta_i ~ N(0, 1)
  
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
##| 29.08030 -2480.21549  4960.43097    30.57148 -2481.11888  4962.23775    30.59223 -2481.13963  4962.27925 

## 3.2 Compute Fully Marginal Likelihood ---------------------------------------

## Step 1: Specify the log_prior function
log_prior <- function(samples_s) {
  a <- as.vector(samples_s[grep("^a\\[\\d+\\]$", names(samples_s))])
  b_1 <- as.vector(samples_s[grep("^b\\[\\d+,\\d+\\]$", names(samples_s))])
  diff <- as.vector(samples_s[grep("^diff\\[\\d+,\\d+\\]$", names(samples_s))])
  
  sum(dlnorm(a, meanlog = 0, sdlog = 1, log = TRUE)) +
    sum(dnorm(b_1, mean = 0, sd = sqrt(1/0.01), log = TRUE)) +
    sum(log(truncdist::dtrunc(diff, "norm", a = 0, b = Inf, mean = 0, sd = sqrt(1/0.01))))
}

## Step 2: Define parameter bounds
lb <- c(rep(0, J), rep(-Inf, J), rep(0, J*(K-2)))
ub <- c(rep(Inf, J), rep(Inf, J), rep(Inf, J*(K-2)))
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
##| Bridge sampling estimate of the log marginal likelihood: -2561.705
##| Estimate obtained in 6 iteration(s) via method "normal".
