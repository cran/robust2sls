% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outlier_tests.R
\name{sumtest}
\alias{sumtest}
\title{Scaling sum proportion test across different cut-offs}
\usage{
sumtest(robust2sls_object, alpha, iteration, one_sided = FALSE)
}
\arguments{
\item{robust2sls_object}{A list of \code{"robust2sls"} objects.}

\item{alpha}{A numeric value between 0 and 1 representing the significance
level of the test.}

\item{iteration}{An integer >= 0 or the character "convergence" that
determines which iteration is used for the test.}

\item{one_sided}{A logical value whether a two-sided test (\code{FALSE})
should be conducted or a one-sided test (\code{TRUE}) that rejects only
when the false outlier detection rate is above its expected value.}
}
\value{
\code{sumtest()} returns a data frame with one row storing the
iteration that was tested, the value of the test statistic (t-test), the
type of the test (one- or two-sided), the corresponding p-value, the
significance level, and whether the null hypothesis is rejected. The data
frame also contains an attribute named \code{"gammas"} that records which
gammas determining the different cut-offs were used in the scaling sum test.
}
\description{
\code{sumtest()} uses the estimations across several cut-offs to test whether
the sum of the deviations between sample and population FODR differ
significantly from its expected value.
\loadmathjax
\mjdeqn{\sum_{k = 1}^{K} \sqrt{n}(\widehat{\gamma}_{c_{k}} - \gamma_{c_{k}}) }{sum k = (1 to K) sqrt(n) (gammahat_{c_{k}} - gamma_{c{k}})}
}
