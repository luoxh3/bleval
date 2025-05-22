
#' @title Information Criteria for Bayesian Latent Variable Models
#'
#' @description
#' Calculates the Deviance Information Criterion (DIC),
#' Watanabe-Akaike Information Criterion (WAIC),
#' and Leave-One-Out Information Criterion (LOOIC)
#' for Bayesian latent variable models based on the log marginal likelihood.
#'
#' @param log_marglik_result A list returned by the `log_marglik()` function, which contains:
#'   \item{log_marglik_point}{A matrix of log marginal likelihood values for each data point.}
#'   \item{log_marglik_postmean}{A vector of log marginal likelihood values computed
#'   using the posterior mean of the parameters.}
#'
#' @param p_dic_version A numeric value indicating which version of the
#'    effective number of parameters (p_D) to use for the DIC calculation.
#'    Version 1 is based on the mean of log marginal likelihood values (default), and
#'    Version 2 is based on the variance of log marginal likelihood values.
#'
#' @return A named vector containing the following information criteria:
#'    \item{p_dic}{The effective number of parameters for DIC.}
#'    \item{elpd_dic}{The expected log predictive density for DIC.}
#'    \item{dic}{The Deviance Information Criterion (DIC), computed as -2 * elpd_dic.}
#'    \item{p_waic}{The effective number of parameters for WAIC.}
#'    \item{elpd_waic}{The expected log predictive density for WAIC.}
#'    \item{waic}{The Watanabe-Akaike Information Criterion (WAIC).}
#'    \item{p_looic}{The effective number of parameters for LOOIC.}
#'    \item{elpd_looic}{The expected log predictive density for LOOIC.}
#'    \item{looic}{The Leave-One-Out Information Criterion (LOOIC).}
#'
#' @importFrom loo loo
#' @importFrom loo waic
#'
#' @export
#'
calc_IC <- function(log_marglik_result, p_dic_version) {

  # Sum of log marginal likelihood values for each replicate
  log_marglik_per_rep <- apply(log_marglik_result$log_marglik_point, 1, sum)

  # Log marginal likelihood based on the posterior mean
  log_marglik_postmean <- sum(log_marglik_result$log_marglik_postmean)

  # Compute p_D (effective number of parameters) for DIC based on the chosen version
  if(p_dic_version == 1) {
    p_dic <-  2 * (log_marglik_postmean - mean(log_marglik_per_rep))
  } else if(p_dic_version == 2){
    p_dic <-  2 * var(log_marglik_per_rep)
  } else{
    # Warn if p_dic_version is invalid and default to version 1
    warning("Invalid p_dic_version: should be either 1 (based on mean) or 2 (based on variance).
            Using version 1 by default.")
    p_dic_version <- 1  # Default to version 1
    p_dic <-  2 * (log_marglik_postmean - mean(log_marglik_per_rep))
  }

  # Deviance Information Criterion (DIC) and expected log predictive density (elpd)
  elpd_dic <- log_marglik_postmean - p_dic

  # Calculate WAIC (Watanabe-Akaike Information Criterion)
  waic_result <- loo::waic(log_marglik_result$log_marglik_point)$estimates
  p_waic <- waic_result[2,1]
  elpd_waic <- waic_result[1,1]
  waic <- waic_result[3,1]

  # Calculate LOOIC (Leave-One-Out Information Criterion)
  loo_result <- loo::loo(log_marglik_result$log_marglik_point)$estimates
  p_looic <- loo_result[2,1]
  elpd_looic <- loo_result[1,1]
  looic <- loo_result[3,1]

  # Return a named vector with all the computed information criteria
  c(p_dic = p_dic, elpd_dic = elpd_dic, dic = -2*elpd_dic,
    p_waic = p_waic, elpd_waic = elpd_waic, waic = waic,
    p_looic = p_looic, elpd_looic = elpd_looic, looic = looic)
}
