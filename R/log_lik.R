
#' @title Log Likelihood per Point and based on Posterior Mean
#'
#' @description
#' This function computes the log likelihood for each data point and
#' the log likelihood based on the posterior mean of the parameters.
#' It is designed for Bayesian latent variable models, where the likelihood
#' is computed by integrating out latent variables using numerical quadrature.
#' The function supports parallel computation to improve efficiency for large datasets.
#'
#'
#' @param samples A matrix or data frame of Bayesian posterior samples of model parameters.
#'    Each row represents a sample, and each column represents a parameter.
#' @param data A list of data, including an element 'N' which indicates the number of units.
#' @param Ngrid Number of grid (quadrature nodes) per dimension.
#' @param lv_mu A list of posterior means for the latent variables.
#'    Each element corresponds to the posterior mean of the latent variables
#'    for a specific unit.
#' @param lv_cov A list of posterior covariance matrix for the latent variables.
#'    Each element corresponds to the posterior covariance matrix of
#'    the latent variables for a specific unit.
#' @param log_joint_i A user-defined function that computes the log joint probability
#'    for a given unit. This function should take the following arguments:
#'    - `samples_s`: A vector of parameter values from a posterior sample.
#'    - `data`: The data list.
#'    - `i`: The index of the unit.
#'    - `Ngrid`: The number of quadrature nodes.
#'    - `nodes`: A matrix of quadrature nodes transformed using the latent variable mean and covariance.
#' @param parallel A logical indicating whether to compute in parallel (default is TRUE).
#' @param n_cores Number of cores to use for parallel computation. Defaults to `detectCores() - 2`.
#' @param packages A character vector of package names to be loaded in the parallel environment.
#'
#' @return A list containing two elements:
#'    \item{loglik_point}{A matrix of log likelihoods for each data point.
#'    Each row corresponds to a posterior sample, and each column corresponds to a unit.}
#'    \item{loglik_postmean}{A vector of log likelihoods computed using the posterior mean
#'    of the parameters. Each element corresponds to a unit.}
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom foreach foreach
#' @importFrom foreach getDoParWorkers
#' @importFrom foreach "%dopar%"
#' @importFrom statmod gauss.quad.prob
#' @importFrom matrixStats logSumExp
#' @importFrom mvtnorm dmvnorm
#' @importFrom extraDistr ddirichlet
#'
#' @export
#'
log_lik <- function(samples, data, Ngrid, lv_mu, lv_cov, log_joint_i,
                    parallel = TRUE, n_cores = detectCores() - 2,
                    packages = c("matrixStats", "statmod", "mvtnorm", "extraDistr")) {

  ## log likelihood for each point ---------------------------------------------

  # Ensure that n_cores is positive and does not exceed the available cores
  if (parallel) {
    if (n_cores <= 0) {
      warning("n_cores must be greater than 0. Defaulting to 1 core.")
      n_cores <- 1
    }

    if (n_cores > detectCores()) {
      warning("The number of cores exceeds available cores. Reducing to available cores.")
      n_cores <- detectCores() - 1
    }

    # Validate that the samples input is either a matrix or a data frame
    if (!is.matrix(samples) && !is.data.frame(samples)) {
      stop("'samples' must be a matrix or data frame.")
    }

    # Validate that the data contains an element 'N'
    if (!"N" %in% names(data)) {
      stop("'data' must contain an element named 'N'.")
    }

    # Start parallel processing using the specified number of cores
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    if (getDoParWorkers() == 0) {
      stop("Parallel backend not registered properly. Please check.")
    }

    loglik_point <- foreach(j = 1:nrow(samples), .combine = rbind,
                            .packages = packages,
                            .export = c("get_quadrature", "log_lik_i",
                                        "data", "Ngrid", "lv_mu", "lv_cov",
                                        "log_joint_i")) %dopar% {
                              sapply(1:data$N, function(i) {
                                bleval::log_lik_i(samples[j, ], data, i, Ngrid,
                                                  lv_mu, lv_cov, log_joint_i)
                              })
                            }
    stopCluster(cl)

  } else {

    # If not using parallel, compute the log likelihood sequentially
    loglik_point <- matrix(nrow = nrow(samples), ncol = data$N)
    for (j in 1:nrow(samples)) {
      for (i in 1:data$N) {
        loglik_point[j, i] <- bleval::log_lik_i(samples[j, ], data, i, Ngrid,
                                                lv_mu, lv_cov, log_joint_i)
      }
    }
  }

  ## log likelihood based on posterior mean ------------------------------------

  samps2_thin_mean <- apply(samples, 2, mean)
  loglik_postmean <- numeric(data$N)
  for (i in 1:data$N) {
    loglik_postmean[i] <- bleval::log_lik_i(samps2_thin_mean, data, i, Ngrid,
                                              lv_mu, lv_cov, log_joint_i)
  }

  # Return a list containing loglik_point and loglik_postmean
  return(list(loglik_point = loglik_point, loglik_postmean = loglik_postmean))
}
