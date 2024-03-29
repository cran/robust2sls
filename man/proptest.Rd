% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outlier_tests.R
\name{proptest}
\alias{proptest}
\title{Proportion test}
\usage{
proptest(robust2sls_object, alpha, iteration, one_sided = FALSE)
}
\arguments{
\item{robust2sls_object}{An object of class \code{"robust2sls"} or a
list of such objects.}

\item{alpha}{A numeric value between 0 and 1 representing the significance
level of the test.}

\item{iteration}{An integer >= 0 or the character "convergence" that
determines which iteration is used for the test.}

\item{one_sided}{A logical value whether a two-sided test (\code{FALSE})
should be conducted or a one-sided test (\code{TRUE}) that rejects only
when the false outlier detection rate is above its expected value.}
}
\value{
\code{proptest()} returns a data frame with the iteration (m) to be
tested, the actual iteration that was tested (generally coincides with the
iteration that was specified to be tested but is the convergent iteration if
the fixed point is tested), the setting of the probability of exceeding the
cut-off (gamma), the type of t-test (one- or two-sided), the value of the
test statistic, its p-value, the significance level \code{alpha}, and the
decision. The number of rows of the data frame corresponds to the length of
the argument \code{robust2sls_object}.
}
\description{
\code{proptest()} conducts a test whether the false outlier detection rate
(FODR) in the sample deviates significantly from its expected value
(population FODR) under the null hypothesis that there are no outliers in the
sample.
}
\details{
See \code{\link[=outlier_detection]{outlier_detection()}} and
\code{\link[=multi_cutoff]{multi_cutoff()}} for creating an object of class
\code{"robust2sls"} or a list thereof.
}
