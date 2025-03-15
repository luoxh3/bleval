
#' @title Log Likelihood for Unit i
#'
#' @param samples_s A vector of Bayesian posterior sample s of model parameters.
#' @param data A list of data, including an element 'N' which indicates the number of units
#' @param i Index of the unit for which to compute the log likelihood.
#' @param Ngrid Number of grid per dimension.
#' @param lv_mu A list of posterior means for the latent variables.
#' @param lv_cov A list of posterior covariance matrix for latent variables.
#' @param log_joint_i A function that computes the log joint probability for unit i.
#'
#' @returns The log likelihood for unit i.
#'
#' @importFrom matrixStats logSumExp
#' @importFrom mvtnorm dmvnorm
#'
#' @export
#'
log_lik_i <- function(samples_s, data, i, Ngrid, lv_mu, lv_cov, log_joint_i) {

  # Extract the posterior mean and covariance matrix for the latent variable for unit i
  lv_mu_i <- lv_mu[[i]]
  lv_cov_i <- lv_cov[[i]]

  # Determine the number of dimensions of the latent variables
  Ndim <- length(lv_mu_i)

  # Compute quadrature nodes and log weights for numerical integration
  nodes_Ndim <- bleval::get_quadrature(Ngrid, Ndim)$nodes
  log_weights_Ndim <- bleval::get_quadrature(Ngrid, Ndim)$log_weights

  # Handle the case when Ndim = 1 (univariate latent variable)
  if (Ndim == 1) {

    # Compute log of the standard normal density for each quadrature point
    # Adjust by the standard deviation (sqrt of variance) for the unit i
    log_std_i <- dnorm(nodes_Ndim, mean = 0, sigma = 1, log = TRUE) - log(sqrt(lv_cov_i))

    # Transform quadrature nodes using the latent variable mean and standard deviation
    nodes <- t(apply(nodes_Ndim, 1, function(z) lv_mu_i + sqrt(lv_cov_i) * z))

  } else {

    # For higher dimensions, compute Cholesky decomposition of the covariance matrix
    lv_Q <- t(chol(lv_cov_i))

    # Compute log of the multivariate normal density for each quadrature point
    # Adjust by the determinant of the Cholesky matrix (covariance scaling)
    log_std_i <- mvtnorm::dmvnorm(nodes_Ndim, mean = rep(0, Ndim),
                                  sigma = diag(Ndim), log = TRUE) - log(det(lv_Q))

    # Transform quadrature nodes using the latent variable mean and the covariance matrix
    nodes <- t(apply(nodes_Ndim, 1, function(z) lv_mu_i + lv_Q %*% z))

  }

  # Compute the log joint probability for unit i using the transformed nodes
  log_joint_i_result <- log_joint_i(samples_s, data, i, Ngrid, nodes)

  # Compute the log likelihood for unit i by summing over the quadrature points
  matrixStats::logSumExp(log_weights_Ndim - log_std_i + log_joint_i_result)

}
