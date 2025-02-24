
#' @title Log Marginal Likelihood for the Model
#'
#' @description
#' Computes the log marginal likelihood of the Bayesian latent variable models
#' using bridge sampling.
#' The function takes posterior samples of the model parameters, a log-likelihood function,
#' and a log-prior function, and then uses the `bridge_sampler` function from the `bridgesampling`
#' package to estimate the log marginal likelihood.
#'
#' @param samples A matrix or data frame of Bayesian posterior samples of model parameters.
#' @param data A list of data, including an element 'N' which indicates the number of persons.
#' @param Ngrid Number of grid per dimension.
#' @param lv_mu A list of posterior means for the latent variables.
#' @param lv_cov A list of posterior covariance matrix for latent variables.
#' @param log_joint_i A function that computes the log joint probability for unit i.
#' @param log_prior A function that computes the log prior probabilities of the model parameters.
#' @param lb named vector with lower bounds for parameters. This is required for
#' the \code{\link[bridgesampling]{bridge_sampler}} function and defines the
#' minimum value for the parameter space over which the marginal likelihood is computed.
#' @param ub named vector with upper bounds for parameters. This is required for
#' the \code{\link[bridgesampling]{bridge_sampler}} function and defines the
#' maximum value for the parameter space over which the marginal likelihood is computed.
#' @param ... Additional arguments passed directly to the `bridge_sampler` function.
#'
#' @return A list containing the result of the bridge sampling, which includes the log marginal
#' likelihood estimate and other information related to the bridge sampling process.
#' The return value is typically an object from the `bridge_sampler` function containing:
#' \itemize{
#'   \item \code{logml}: estimate of log marginal likelihood.
#'   \item \code{niter}: number of iterations of the iterative updating scheme.
#' }
#'
#' @importFrom statmod gauss.quad.prob
#' @importFrom matrixStats logSumExp
#' @importFrom mvtnorm dmvnorm
#' @importFrom bridgesampling bridge_sampler
#'
#' @export
#'
#' @examples
#' \donttest{
#' log_marg_lik(...)  # need another file
#' }
#'
log_marg_lik <- function(samples, data, Ngrid, lv_mu, lv_cov, log_joint_i,
                         log_prior, lb, ub, ...) {

  # Define log_posterior function for bridge_sampler ---------------------------
  # Note: use closure to capture lv_mu, lv_cov, Ngrid and log_joint_i
  log_posterior <- function(samples_s, data) {

    N <- data$N
    loglik_i_all <- numeric(N)
    for (i in 1:N) {
      # Compute log likelihood for each individual
      # passing lv_mu, lv_cov and Ngrid explicitly
      loglik_i_all[i] <- blvmeval::log_lik_i(samples_s, data, i, Ngrid,
                                             lv_mu, lv_cov, log_joint_i)
    }
    # Return the sum of the log likelihoods and the log prior
    sum(loglik_i_all) + log_prior(samples_s)
  }

  # Use bridge_sampler to compute log marginal likelihood ----------------------
  bridge_result <- bridgesampling::bridge_sampler(samples = samples,
                                                  log_posterior = log_posterior,
                                                  data = data,
                                                  lb = lb, ub = ub, ...)
  # Return the result of the bridge sampling
  return(bridge_result)

}
