% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{case_resampling}
\alias{case_resampling}
\title{Uses nonparametric case resampling for standard errors of parameters and
gauge}
\usage{
case_resampling(robust2sls_object, R, coef = NULL, m = NULL, parallel = FALSE)
}
\arguments{
\item{robust2sls_object}{An object of class \code{"robust2sls"}.}

\item{R}{An integer specifying the number of resamples.}

\item{coef}{A numeric or character vector specifying which structural
coefficient estimates should be recorded across bootstrap replications.
\code{NULL} means all coefficients are recorded.}

\item{m}{A single numeric or vector of integers specifying for which
iterations the bootstrap statistics should be calculated. \code{NULL} means
they are calculated for all iterations that were also done in the original
robust2sls_object. Character \code{"convergence"} means all bootstrap samples
are run until they converge and the statistics of the first convergent
iteration is recorded.}

\item{parallel}{A logical value indicating whether to run the bootstrap
sampling in parallel or sequentially. See Details.}
}
\value{
\code{case_resampling} returns an object of class
\code{"r2sls_boot"}. This is a list with three named elements. \code{$boot}
stores the bootstrap results as a data frame. The columns record the
different test statistics, the iteration \code{m}, and the number of the
resample, \code{r}. The values corresponding to the original data is stored
as \code{r = 0}. \code{$resamples} is a list of length \code{R} that stores
the indices for each specific resample. \code{$original} stores the original
\code{robust2sls_object} based on which the bootstrapping was done.
}
\description{
Uses nonparametric case resampling for standard errors of parameters and
gauge
}
\details{
Argument \code{parallel} allows for parallel computing using the
\link[foreach]{foreach} package, so the user has to register a parallel
backend before invoking this command.

Argument \code{coef} is useful if the model includes many controls whose
parameters are not of interest. This can reduce the memory space needed to
store the bootstrap results.
}
