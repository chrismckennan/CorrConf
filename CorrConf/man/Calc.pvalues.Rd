% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EstimateC_complete.R
\name{Calc.pvalues}
\alias{Calc.pvalues}
\title{Estimate beta and compute p-values under the model Y ~ MN\{ beta \%*\% X' + Gamma \%*\% Z' + L \%*\% C', diag(delta_1^2,...,delta_p^2), V(tau) \}}
\usage{
Calc.pvalues(
  Y,
  B = NULL,
  X,
  Z = NULL,
  C = NULL,
  tau = NULL,
  A.ine = NULL,
  D.ker = NULL
)
}
\arguments{
\item{Y}{A \code{p} x \code{n} data matrix (# of genes x # of samples)}

\item{B}{A list of \code{n} x \code{n} positive semi-definite matrices. Default is NULL}

\item{X}{The \code{n} x \code{d} covariates of interest}

\item{Z}{A \code{n} x \code{r} matrix of nuisance covariates. Default is NULL}

\item{C}{The estimated \code{n} x \code{K} matrix of latent factors (see EstimateC_complete). Default is NULL}

\item{tau}{The starting point for tau. Default is NULL and the user need not supply one}

\item{A.ine}{A #inequality constaints x b matrix, where b = length(B) = num. variance components. \code{A.ine \%*\% tau >= 0}, where tau = (tau_1^2,...,tau_b^2) are the variance multipliers. The default is all variance multipliers are greater than 0.}

\item{D.ker}{No need to specify this}
}
\value{
A list \item{t}{The \code{p} x \code{d} matrix of t-scores} \item{p}{The \code{p} x \code{d} matrix of p-values} \item{Beta.hat}{The estimate for beta, a \code{p} x \code{d} matrix} \item{Delta.hat}{A \code{p}-vector of residual variances} \item{tau}{The estimate for tau}
}
\description{
Estimate beta and compute p-values under the model Y ~ MN\{ beta \%*\% X' + Gamma \%*\% Z' + L \%*\% C', diag(delta_1^2,...,delta_p^2), V(tau) \}
}
