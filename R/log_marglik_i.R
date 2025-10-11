
#' @title Log Marginal Likelihood for Unit i
#'
#' @description
#' This function computes the log marginal likelihood for a specific unit by
#' integrating out the latent variables using numerical quadrature.
#' It is a helper function used internally by `log_marglik()` to compute the
#' log marginal likelihood for each unit.
#'
#' @param samples_s A vector of Bayesian posterior samples of model parameters.
#' @param data A list of data, including an element 'N' which indicates the number of units
#' @param i Index of the unit for which to compute the log marginal likelihood.
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
#' @param ... Additional arguments to the log_joint_i function.
#'
#' @returns The log marginal likelihood for unit i.
#'
#' @importFrom matrixStats logSumExp
#' @importFrom mvtnorm dmvnorm
#'
#' @export
#'
log_marglik_i <- function(samples_s, data, i, Ngrid, lv_mu, lv_cov, log_joint_i, ...) {

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
    lv_mu_i <- as.numeric(lv_mu_i)
    lv_cov_i <- as.numeric(lv_cov_i)

    # Compute log of the standard normal density for each quadrature point
    # Adjust by the standard deviation (sqrt of variance) for the unit i
    log_std_i <- t( dnorm(nodes_Ndim, mean = 0, sd = 1, log = TRUE) - log(sqrt(lv_cov_i)) )

    # Transform quadrature nodes using the latent variable mean and standard deviation
    nodes <- lv_mu_i + sqrt(lv_cov_i) * nodes_Ndim

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

  # Compute the log joint density for unit i using the transformed nodes
  log_joint_i_result <- log_joint_i(samples_s, data, i, Ngrid, nodes, ...)

  # Compute the log marginal likelihood for unit i by summing over the quadrature points
  matrixStats::logSumExp(log_weights_Ndim - log_std_i + log_joint_i_result)

}
