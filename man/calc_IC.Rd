% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_IC.R
\name{calc_IC}
\alias{calc_IC}
\title{Information Criteria for Bayesian Latent Variable Models}
\usage{
calc_IC(log_marglik_result, p_dic_version)
}
\arguments{
\item{log_marglik_result}{A list returned by the \code{log_marglik()} function, which contains:
\item{log_marglik_point}{A matrix of log marginal likelihood values for each data point.}
\item{log_marglik_postmean}{A vector of log marginal likelihood values computed
using the posterior mean of the parameters.}}

\item{p_dic_version}{A numeric value indicating which version of the
effective number of parameters (p_D) to use for the DIC calculation.
Version 1 is based on the mean of log marginal likelihood values (default), and
Version 2 is based on the variance of log marginal likelihood values.}
}
\value{
A named vector containing the following information criteria:
\item{p_dic}{The effective number of parameters for DIC.}
\item{elpd_dic}{The expected log predictive density for DIC.}
\item{dic}{The Deviance Information Criterion (DIC), computed as -2 * elpd_dic.}
\item{p_waic}{The effective number of parameters for WAIC.}
\item{elpd_waic}{The expected log predictive density for WAIC.}
\item{waic}{The Watanabe-Akaike Information Criterion (WAIC).}
\item{p_looic}{The effective number of parameters for LOOIC.}
\item{elpd_looic}{The expected log predictive density for LOOIC.}
\item{looic}{The Leave-One-Out Information Criterion (LOOIC).}
}
\description{
Calculates the Deviance Information Criterion (DIC),
Watanabe-Akaike Information Criterion (WAIC),
and Leave-One-Out Information Criterion (LOOIC)
for Bayesian latent variable models based on the log marginal likelihood.
}
