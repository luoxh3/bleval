
#' @title Log Fully Marginal Likelihood for the Model
#'
#' @description
#' This function computes the log fully marginal likelihood of the Bayesian
#' latent variable models using bridge sampling.
#' The log fully marginal likelihood is computed by integrating out both
#' the latent variables and the model parameters.
#' The function uses the `bridge_sampler` function from the `bridgesampling` package
#' to estimate the log fully marginal likelihood.
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
#' @param log_prior A user-defined function that computes the log prior density
#'    of the model parameters. This function should take the following argument:
#'    - `samples_s`: A vector of parameter values from a posterior sample.
#' @param lb A named vector with lower bounds for parameters. This is required for
#'    the \code{\link[bridgesampling]{bridge_sampler}} function and defines the
#'    minimum value for the parameter space over which the fully marginal likelihood is computed.
#' @param ub A named vector with upper bounds for parameters. This is required for
#'    the \code{\link[bridgesampling]{bridge_sampler}} function and defines the
#'    maximum value for the parameter space over which the fully marginal likelihood is computed.
#' @param ... Additional arguments passed directly to the `bridge_sampler` function.
#'
#' @return A list containing the result of the bridge sampling, which includes the log fully marginal
#'    likelihood estimate and other information related to the bridge sampling process.
#'    The return value is typically an object from the `bridge_sampler` function containing:
#'    \itemize{
#'      \item \code{logfml}: estimate of log fully marginal likelihood.
#'      \item \code{niter}: number of iterations of the iterative updating scheme.
#' }
#'
#' @importFrom statmod gauss.quad.prob
#' @importFrom matrixStats logSumExp
#' @importFrom mvtnorm dmvnorm
#' @importFrom bridgesampling bridge_sampler
#'
#' @export
#'
log_fmarglik <- function(samples, data, Ngrid, lv_mu, lv_cov, log_joint_i,
                         log_prior, lb, ub, ...) {

  # Define log_posterior function for bridge_sampler ---------------------------
  # Note: use closure to capture lv_mu, lv_cov, Ngrid and log_joint_i
  log_posterior <- function(samples_s, data) {

    N <- data$N
    log_marglik_i_all <- numeric(N)
    for (i in 1:N) {
      # Compute log marginal likelihood for each unit
      # passing lv_mu, lv_cov and Ngrid explicitly
      log_marglik_i_all[i] <- bleval::log_marglik_i(samples_s, data, i, Ngrid,
                                                    lv_mu, lv_cov, log_joint_i)
    }
    # Return the sum of the log marginal likelihoods and the log prior density
    sum(log_marglik_i_all) + log_prior(samples_s)
  }

  # Use bridge_sampler to compute log fully marginal likelihood ----------------
  bridge_result <- bridgesampling::bridge_sampler(samples = samples,
                                                  log_posterior = log_posterior,
                                                  data = data,
                                                  lb = lb, ub = ub, ...)
  # Return the result of the bridge sampling
  return(bridge_result)

}

