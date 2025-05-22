
#' @title Standard Quadrature Nodes and Weights
#'
#' @description
#' Generates the quadrature nodes and weights using Gaussian quadrature
#' for a specified number of grid points and dimensions.
#' The quadrature is based on a standard normal distribution (mean = 0, sd = 1).
#' The function is used internally by `log_marglik_i()` to perform
#' numerical integration over latent variables.
#'
#' @param Ngrid Number of grid points (quadrature nodes) per dimension.
#' @param Ndim Number of dimensions for the latent variables.
#'    This corresponds to the number of latent variables in the model.
#'
#' @return A list containing the following elements:
#'    \item{nodes}{A matrix of quadrature nodes with dimensions
#'    \eqn{N_{\text{nodes}} \times N_{\text{dim}}}{Nnodes × Ndim}, where
#'    \eqn{N_{\text{nodes}} = \text{Ngrid}^{\text{Ndim}}}{Nnodes = Ngrid^Ndim}.
#'    Each column represents a specific combination of nodes from the 1D Gaussian
#'    quadrature expanded into the multidimensional grid using `expand.grid`.}
#'    \item{weights}{A vector of quadrature weights corresponding to the nodes.}
#'    \item{log_weights}{A vector of the log of the weights,
#'    typically used in logarithmic space computations.}
#'
#' @importFrom statmod gauss.quad.prob
#'
#' @export
#'
#' @examples
#' get_quadrature(5, 2)  # Generate quadrature nodes and weights for 2D
#'                       # with 5 grid points per dimension.
#'
get_quadrature <- function(Ngrid, Ndim) {

  # Generate standard Gaussian quadrature for 1D
  std_quad <- statmod::gauss.quad.prob(Ngrid, "normal", mu = 0, sigma = 1)
  nodes_1d <- std_quad$nodes
  weights_1d <- std_quad$weights

  # Create Ndim-dimensional grid by expanding the 1D nodes
  nodes_nd <- as.matrix(do.call(expand.grid, rep(list(nodes_1d), Ndim)))

  # Compute the Ndim-dimensional weights as the product of individual 1D weights
  weights_nd_temp <- as.matrix(do.call(expand.grid, rep(list(weights_1d), Ndim)))
  weights_nd <- apply(weights_nd_temp, 1, prod)
  log_weights_nd <- log(weights_nd)

  # Return the quadrature nodes, weights, and log_weights
  list(nodes = nodes_nd, weights = weights_nd, log_weights = log_weights_nd)
}
