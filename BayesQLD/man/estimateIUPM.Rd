% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfunctions.R
\name{estimateIUPM}
\alias{estimateIUPM}
\title{Estimate IUPM from data}
\usage{
estimateIUPM(
  data,
  well.state.colID,
  dilution.colID,
  nBurnin = 1000,
  nSampling = 5000,
  nThin = 2,
  nChains = 4,
  probability = 0.95
)
}
\arguments{
\item{data}{data frame with data}

\item{well.state.colID}{name or number of the column containing well state data}

\item{dilution.colID}{name or number of the column with dilution factors}

\item{nBurnin}{number of burn-in iterations (default is 1000, the recommended value)}

\item{nSampling}{number of sampling iterations per chain (default is 12500, the recommended value)}

\item{nThin}{chain thinning value (default is 5, the recommended value)}

\item{nChains}{number of chains (default is 4, the recommended value)}

\item{probability}{the probability of the outer HPD quantile to report (default is 0.95)}
}
\value{
a list containing the highest-posterior density quantiles (the outer, inter-quartile range, and mode) of IUPM estimates and MCMC convergence diagnostics
}
\description{
Fits a Bayesian binomial model to data, estimating IUPM. The data should be in a data frame that contains at least two columns: (1) positive or negative status and (2) the dilution factor for each well. Positive wells must be 1, negative 0. Each row of the data frame is a well. The dilution factor is the fraction of the initial concentration. If the initial number of cells is 1 million, the result is IUPM. If the initial number is different, however, the reported values must be corrected by dividing the quantiles of the reported value by the initial number of cells.
}
\references{
\insertRef{fazekas82a}{BayesQLD}
}
