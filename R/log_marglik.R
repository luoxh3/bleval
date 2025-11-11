
#' @title Log Marginal Likelihood per Point and based on Posterior Mean
#'
#' @description
#' This function computes the log marginal likelihood for each data point and
#' the log marginal likelihood based on the posterior mean of the parameters.
#' It is designed for Bayesian latent variable models, where the marginal likelihood
#' is computed by integrating out latent variables using numerical quadrature.
#' The function supports parallel computation to improve efficiency for large datasets.
#'
#'
#' @param samples A matrix or data frame of Bayesian posterior samples of model parameters.
#'    Each row represents a sample, and each column represents a parameter.
#' @param data A list of data, including an element 'N' which indicates the number of units.
#' @param Ngrid Number of grid points (quadrature nodes) per dimension.
#' @param lv_mu A list of posterior means for the latent variables.
#'    Each element corresponds to the posterior mean of the latent variables
#'    for a specific unit.
#' @param lv_cov A list of posterior covariance matrix for the latent variables.
#'    Each element corresponds to the posterior covariance matrix of
#'    the latent variables for a specific unit.
#' @param log_joint_i A user-defined function that computes the log joint density
#'    for a given unit. This function should take the following arguments:
#'    - `samples_s`: A vector of parameter values from a posterior sample.
#'    - `data`: The data list.
#'    - `i`: The index of the unit.
#'    - `Ngrid`: Number of grid points (quadrature nodes) per dimension.
#'    - `nodes`: A matrix of quadrature nodes transformed using the latent variable mean and covariance.
#' @param parallel A logical indicating whether to compute in parallel (default is TRUE).
#' @param n_cores Number of cores to use for parallel computation. Defaults to `detectCores() - 2`.
#' @param packages A character vector of additional package names to be loaded in the parallel environment.
#' @param exports A character vector of additional function/object names to export to parallel workers.
#' @param ... Additional arguments to the log_joint_i function.
#'
#' @return A list containing two elements:
#'    \item{log_marglik_point}{A matrix of log marginal likelihoods for each data point.
#'    Each row corresponds to a posterior sample, and each column corresponds to a unit.}
#'    \item{log_marglik_postmean}{A vector of log marginal likelihoods computed using the posterior mean
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
#' @importFrom truncdist dtrunc
#'
#' @export
#'
log_marglik <- function(samples, data, Ngrid, lv_mu, lv_cov, log_joint_i,
                    parallel = TRUE, n_cores = detectCores() - 2,
                    packages = NULL, exports = NULL, ...) {

  ## log marginal likelihood for each point ------------------------------------
  required_packages <- c("matrixStats", "statmod", "mvtnorm", "extraDistr", "truncdist")
  if (!is.null(packages)) {
    required_packages <- unique(c(required_packages, packages))
  }
  required_exports <- c("get_quadrature", "log_marglik_i",
                        "data", "Ngrid", "lv_mu", "lv_cov",
                        "log_joint_i")
  if (!is.null(exports)) {
    required_exports <- unique(c(required_exports, exports))
  }
  
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
      stop("'data' must contain an element named 'N' which indicates the number of units.")
    }

    if (.Platform$OS.type == "windows"){
      # Start parallel processing using the specified number of cores
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)

      if (getDoParWorkers() == 0) {
        stop("Parallel backend not registered properly. Please check.")
      }

      log_marglik_point <- foreach(j = 1:nrow(samples), .combine = rbind,
                                  .packages = required_packages,
                                  .export = required_exports) %dopar% {
                                    sapply(1:data$N, function(i) {
                                      bleval::log_marglik_i(samples[j, ], data, i, Ngrid,
                                                            lv_mu, lv_cov, log_joint_i, ...)
                                    })
                                  }
      stopCluster(cl)
    } else {
      missing_exports <- setdiff(required_exports, ls(envir = .GlobalEnv))
      if (length(missing_exports) > 0) {
        warning("The following functions/objects are not found in global environment and may cause errors: ",
                paste(missing_exports, collapse = ", "))
      }
      
      log_marglik_list <- mclapply(1:nrow(samples), function(j) {
        sapply(1:data$N, function(i) {
          bleval::log_marglik_i(samples[j, ], data, i, Ngrid, lv_mu, lv_cov, log_joint_i, ...)
        })
      }, mc.cores = n_cores)
      log_marglik_point <- do.call(rbind, log_marglik_list)
    }

  } else {

    # If not using parallel, compute the log marginal likelihood sequentially
    log_marglik_point <- matrix(nrow = nrow(samples), ncol = data$N)
    for (j in 1:nrow(samples)) {
      for (i in 1:data$N) {
        log_marglik_point[j, i] <- bleval::log_marglik_i(samples[j, ], data, i, Ngrid,
                                                         lv_mu, lv_cov, log_joint_i, ...)
      }
    }
  }

  ## log marginal likelihood based on posterior mean ---------------------------

  samps2_thin_mean <- apply(samples, 2, mean)
  log_marglik_postmean <- numeric(data$N)
  for (i in 1:data$N) {
    log_marglik_postmean[i] <- bleval::log_marglik_i(samps2_thin_mean, data, i, Ngrid,
                                                     lv_mu, lv_cov, log_joint_i, ...)
  }

  # Return a list containing log_marglik_point and log_marglik_postmean
  return(list(log_marglik_point = log_marglik_point,
              log_marglik_postmean = log_marglik_postmean))
}
