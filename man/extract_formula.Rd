% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{extract_formula}
\alias{extract_formula}
\title{Extract the elements of ivreg formula}
\usage{
extract_formula(formula)
}
\arguments{
\item{formula}{A formula for the \code{\link[ivreg]{ivreg}} function, i.e. in
format \code{y ~ x1 + x2 | z1 + z2}.}
}
\value{
\code{extract_formula} returns a list with five named components,
each of which is a character vector: \code{$y_var} refers to the dependent
variable, \code{$x1_var} to the exogenous regressors, \code{$x2_var} to the
endogenous regressors, \code{$z1_var} to the exogenous regressors (which have
to be included again as instruments and hence coincide with \code{$x1_var}),
and \code{$z2_var} refers to the outside instruments.
}
\description{
\code{extract_formula} takes a formula object for \code{\link[ivreg]{ivreg}},
i.e. in a format of \code{y ~ x1 + x2 | x1 + z2} and extracts the different
elements in a list. Each element is a character vector storing the different
types of regressors. Element \code{y_var} refers to the dependent variable,
\code{x1_var} to the exogenous regressors, \code{x2_var} to the endogenous
regressors, \code{z1_var} to the exogenous regressors (which have to be
included again as instruments and hence coincide with \code{x1_var}), and
\code{z2_var} refers to the outside instruments.
}
\keyword{internal}
