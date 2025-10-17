## This R script demonstrates Bayesian evaluation of a Generalized Partial Credit Model (GPCM).
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

# Specify the data-generating GPCM:
a <- c(1.0, 1.0, 1.5, 1.5, 2.0, 2.0) # Item discrimination parameters
b <- matrix(rep(c(-1.0, -0.5, 0.5, 1.5), times = J), 
            nrow = J, ncol = K-1, byrow = TRUE) + 
  runif(J, -0.5, 0.5) # Item difficulty parameters
b_mat <- d_mat <- matrix(0, nrow = J, ncol = K)
b_mat[,2:K] <- b*a

for (j in 1:J)
  d_mat[j, ] <- -cumsum(b_mat[j,])

# Generate GPCM dataset
data <- simdata(a = a, d = d_mat, N = N, itemtype = 'gpcm')

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
      for (k in 1:K) {
        comp[i,j,k] <- a[j] * ((k-1) * theta[i] - sum_b[j,k])
        exp_comp[i,j,k] <- exp(comp[i,j,k])
      }
      sum_comp[i,j] <- sum(exp_comp[i,j,1:K])

      
      for (k in 1:K) {
        p[i,j,k] <- exp_comp[i,j,k] / sum_comp[i,j]
      }
      
      y[i,j] ~ dcat(p[i,j,1:K])
    }
  }
  
  for (j in 1:J) {
    sum_b[j,1] <- 0
    for (k in 2:K) {
      sum_b[j,k] <- sum(b[j,1:(k-1)])
    }
  }
  
  # Priors -------------------------------------
  for (j in 1:J) {
    a[j] ~ dlnorm(0, 1)
    for (k in 1:(K-1)) {
      b[j,k] ~ dnorm(0, 0.01)
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
                     variable.names = c("a", "b", "theta"),
                     n.iter = 25000, thin = 10)

# Combine posterior samples from all chains into a single matrix
samps <- do.call(rbind, post)

# Extract posterior samples for model parameters only (excluding latent variables)
# 30 model parameters
pars_vector <- paste("a[",1:J,"]", sep="")
for (j in 1:J)
  for (k in 1:(K-1))
    pars_vector <- c(pars_vector, paste("b[",j,",",k,"]", sep=""))

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
  b <- matrix(samples_s[grep("^b\\[\\d+,\\d+\\]$", names(samples_s))],
              nrow = J, ncol = K-1, byrow = TRUE)
  
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
  # In GPCM, category probabilities are computed via softmax transformation 
  # of the linear components: a_j * ((k-1) * theta_i + d_jk)
  
  prob <- matrix(NA, nrow = Ngrid^1, ncol = K)
  log_con_i_j <- matrix(NA, nrow = Ngrid^1, ncol = J)
  
  for (j in 1:J) {
    # Compute category probabilities using GPCM formulation
    d <- cumsum(c(0, -b[j,]))
    d_mat <- matrix(rep(d, times = Ngrid^1),
                    nrow = Ngrid^1, ncol = K, byrow = TRUE)
    ak_mat <- matrix(rep(0:(K-1), times = Ngrid^1),
                     nrow = Ngrid^1, ncol = K, byrow = TRUE)
    numerator <- exp(a[j] * (as.numeric(nodes) * ak_mat + d_mat))
    denominator <- rowSums(numerator)
    prob <- numerator/denominator
    
    # Compute log-density for each observation at each quadrature point
    log_con_i_j[,j] <- extraDistr::dcat(y_i[j], prob, log = TRUE)
  }
  
  # Sum over observations within the unit (conditional independence)
  log_con_i <- rowSums(log_con_i_j)
  
  # ==================================================
  # STEP 4: Compute log density for latent variables
  # ==================================================
  # This corresponds to: log p(η_j | θ)
  
  # In GPCM, we have a normal distribution for person i's latent trait: theta_i ~ N(0, 1)
  
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
##    p_dic    elpd_dic         dic      p_waic   elpd_waic        waic     p_looic  elpd_looic       looic
## 30.44384 -2407.10130  4814.20260    30.07185 -2407.09580  4814.19160    30.09403 -2407.11798  4814.23596

## 3.2 Compute Fully Marginal Likelihood ---------------------------------------

## Step 1: Specify the log_prior function
log_prior <- function(samples_s) {
  a <- as.vector(samples_s[grep("^a\\[\\d+\\]$", names(samples_s))])
  b <- matrix(samples_s[grep("^b\\[\\d+,\\d+\\]$", names(samples_s))],
              nrow = J, ncol = K-1, byrow = TRUE)
  
  sum(dlnorm(a, meanlog = 0, sdlog = 1, log = TRUE)) 
  sum(dnorm(b, mean = 0, sd = sqrt(1/0.01), log = TRUE))
}

## Step 2: Define parameter bounds
lb <- c(rep(0, J), rep(-Inf, J * (K-1)))
ub <- c(rep(Inf, J), rep(Inf, J * (K-1)))
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
##| Bridge sampling estimate of the log marginal likelihood: -2479.794
##| Estimate obtained in 8 iteration(s) via method "normal".
