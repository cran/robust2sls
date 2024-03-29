% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{nonparametric_resampling}
\alias{nonparametric_resampling}
\title{Nonparametric resampling from a data frame}
\usage{
nonparametric_resampling(df, resample)
}
\arguments{
\item{df}{Data frame containing observations to be sampled from.}

\item{resample}{A vector of indices that extract the observations from the
data frame.}
}
\value{
\code{nonparametric_resampling} returns a data frame containing the
observations of the resample.
}
\description{
Nonparametric resampling from a data frame
}
\details{
The input to the \code{resample} argument could for example be generated as
one of the elements in the list generated by the command
\link{nonparametric}.

The input to the \code{df} argument would be the original data frame for case
resampling. For error/residual resampling, it would be a data frame
containing the residuals from the model.
}
