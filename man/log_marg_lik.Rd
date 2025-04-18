% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_marg_lik.R
\name{log_marg_lik}
\alias{log_marg_lik}
\title{Log Marginal Likelihood for the Model}
\usage{
log_marg_lik(
  samples,
  data,
  Ngrid,
  lv_mu,
  lv_cov,
  log_joint_i,
  log_prior,
  lb,
  ub,
  ...
)
}
\arguments{
\item{samples}{A matrix or data frame of Bayesian posterior samples of model parameters.
Each row represents a sample, and each column represents a parameter.}

\item{data}{A list of data, including an element 'N' which indicates the number of units.}

\item{Ngrid}{Number of grid points (quadrature nodes) per dimension.}

\item{lv_mu}{A list of posterior means for the latent variables.
Each element corresponds to the posterior mean of the latent variables
for a specific unit.}

\item{lv_cov}{A list of posterior covariance matrix for the latent variables.
Each element corresponds to the posterior covariance matrix of
the latent variables for a specific unit.}

\item{log_joint_i}{A user-defined function that computes the log joint probability
for a given unit. This function should take the following arguments:
\itemize{
\item \code{samples_s}: A vector of parameter values from a posterior sample.
\item \code{data}: The data list.
\item \code{i}: The index of the unit.
\item \code{Ngrid}: Number of grid points (quadrature nodes) per dimension.
\item \code{nodes}: A matrix of quadrature nodes transformed using the latent variable mean and covariance.
}}

\item{log_prior}{A user-defined function that computes the log prior probabilities
of the model parameters. This function should take the following argument:
\itemize{
\item \code{samples_s}: A vector of parameter values from a posterior sample.
}}

\item{lb}{A named vector with lower bounds for parameters. This is required for
the \code{\link[bridgesampling]{bridge_sampler}} function and defines the
minimum value for the parameter space over which the marginal likelihood is computed.}

\item{ub}{A named vector with upper bounds for parameters. This is required for
the \code{\link[bridgesampling]{bridge_sampler}} function and defines the
maximum value for the parameter space over which the marginal likelihood is computed.}

\item{...}{Additional arguments passed directly to the \code{bridge_sampler} function.}
}
\value{
A list containing the result of the bridge sampling, which includes the log marginal
likelihood estimate and other information related to the bridge sampling process.
The return value is typically an object from the \code{bridge_sampler} function containing:
\itemize{
\item \code{logml}: estimate of log marginal likelihood.
\item \code{niter}: number of iterations of the iterative updating scheme.
}
}
\description{
This function computes the log marginal likelihood of the Bayesian
latent variable models using bridge sampling.
The log marginal likelihood is computed by integrating out both
the latent variables and the model parameters.
The function uses the \code{bridge_sampler} function from the \code{bridgesampling} package
to estimate the log marginal likelihood.
}
