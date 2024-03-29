% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.mrhevo.R
\name{weighted.median.boot}
\alias{weighted.median.boot}
\title{Bootstrap standard error for weighted median estimate of ratio of gamma to
alpha, where gamma and alpha are independent Gaussian variates with
estimates gamma_hat, alpha_hat and respective standard errors se.gamma_hat,
se.alpha_hat.}
\usage{
weighted.median.boot(alpha_hat, gamma_hat, se.alpha_hat, se.gamma_hat, weights)
}
\arguments{
\item{alpha_hat}{Estimate of Gaussian alpha.}

\item{gamma_hat}{Estimate of Gaussian gamma.}

\item{se.alpha_hat}{Standard error of Gaussian alpha.}

\item{se.gamma_hat}{Standard error of Gaussian gamma.}
}
\value{
Standard error of weighted median estimate.
}
\description{
Bootstrap standard error for weighted median estimate of ratio of gamma to
alpha, where gamma and alpha are independent Gaussian variates with
estimates gamma_hat, alpha_hat and respective standard errors se.gamma_hat,
se.alpha_hat.
}
