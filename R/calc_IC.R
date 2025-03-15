
#' @title Information Criteria for Bayesian Latent Variable Models
#'
#' @description
#' Calculates the Deviance Information Criterion (DIC),
#' Widely Applicable Information Criterion (WAIC),
#' and Leave-One-Out #' Information Criterion (LOOIC)
#' for Bayesian latent variable models based on the log likelihood.
#'
#' @param log_lik_result A list returned by the `log_lik()` function, which contains:
#'   \item{loglik_point}{A matrix of log-likelihood values for each data point.}
#'   \item{loglik_postmean}{A vector of log-likelihood values computed using the posterior mean of the parameters.}
#'
#' @param p_dic_version A numeric value indicating which version of the
#'    effective number of parameters (p_D) to use for the DIC calculation.
#'    Version 1 is based on the mean of log-likelihood values (default), and
#'    Version 2 is based on the variance of log-likelihood values.
#'
#' @return A named vector containing the following information criteria:
#'    \item{p_dic}{The effective number of parameters for DIC.}
#'    \item{elpd_dic}{The expected log pointwise predictive density for DIC.}
#'    \item{dic}{The Deviance Information Criterion (DIC), computed as -2 * elpd_dic.}
#'    \item{p_waic}{The effective number of parameters for WAIC.}
#'    \item{elpd_waic}{The expected log pointwise predictive density for WAIC.}
#'    \item{waic}{The Widely Applicable Information Criterion (WAIC).}
#'    \item{p_looic}{The effective number of parameters for LOOIC.}
#'    \item{elpd_looic}{The expected log pointwise predictive density for LOOIC.}
#'    \item{looic}{The Leave-One-Out Information Criterion (LOOIC).}
#'
#' @importFrom loo loo
#' @importFrom loo waic
#'
#' @export
#'
calc_IC <- function(log_lik_result, p_dic_version) {

  # Sum of log-likelihood values for each replicate
  loglik_per_rep <- apply(log_lik_result$loglik_point, 1, sum)

  # Log-likelihood based on the posterior mean
  loglik_postmean <- sum(log_lik_result$loglik_postmean)

  # Compute p_D (effective number of parameters) for DIC based on the chosen version
  if(p_dic_version == 1) {
    p_dic <-  2 * (loglik_postmean - mean(loglik_per_rep))
  } else if(p_dic_version == 2){
    p_dic <-  2 * var(loglik_per_rep)
  } else{
    # Warn if p_dic_version is invalid and default to version 1
    warning("Invalid p_dic_version: should be either 1 (based on mean) or 2 (based on variance).
            Using version 1 by default.")
    p_dic_version <- 1  # Default to version 1
    p_dic <-  2 * (loglik_postmean - mean(loglik_per_rep))
  }

  # Deviance Information Criterion (DIC) and expected log predictive density (elpd)
  elpd_dic <- loglik_postmean - p_dic

  # Calculate WAIC (Widely Applicable Information Criterion)
  waic_result <- loo::waic(log_lik_result$loglik_point)$estimates
  p_waic <- waic_result[2,1]
  elpd_waic <- waic_result[1,1]
  waic <- waic_result[3,1]

  # Calculate LOOIC (Leave-One-Out Information Criterion)
  loo_result <- loo::loo(log_lik_result$loglik_point)$estimates
  p_looic <- loo_result[2,1]
  elpd_looic <- loo_result[1,1]
  looic <- loo_result[3,1]

  # Return a named vector with all the computed information criteria
  c(p_dic = p_dic, elpd_dic = elpd_dic, dic = -2*elpd_dic,
    p_waic = p_waic, elpd_waic = elpd_waic, waic = waic,
    p_looic = p_looic, elpd_looic = elpd_looic, looic = looic)
}
