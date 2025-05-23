
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R Package “robust2sls”

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/jkurle/robust2sls/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jkurle/robust2sls?branch=master)
[![](https://www.r-pkg.org/badges/version/robust2sls?color=blue)](https://cran.r-project.org/package=robust2sls)
[![R-CMD-check](https://github.com/jkurle/robust2sls/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jkurle/robust2sls/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of robust2sls is to provide easy-to-use tools for
outlier-robust inference and outlier testing in two-stage least squares
(2SLS) models.

## Installation

You can install the released version from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("robust2sls")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jkurle/robust2sls")
```

## Introduction

For a detailed introduction to the model framework, the different
trimmed 2SLS algorithms, and examples, see the vignette *Introduction to
the robust2sls Package*.

``` r
utils::vignette("overview", package = "robust2sls")
```

## Note about Versions of Dependencies

The `Depends:` and `Suggests:` fields in the `DESCRIPTION` file have no
minimum or maximum version because I cannot test how far the package is
compatible with older versions of the dependencies. However, I also did
not want to require the versions that were used during the development
to not force users to update their packages and potentially break their
other existing code.

The following table lists the version of each package that was used in
the development of the *robust2sls* package.

For your information, *robust2sls* was developed under the following
versions:

- R: version 4.1.1 (2021-08-10)
- Imports:
  - AER: v1.2-9
  - doRNG: v1.8.2
  - foreach: v1.5.1
  - pracma: v2.3.3
- Suggests:
  - doFuture: v0.12.0
  - doParallel: v1.0.16
  - future: v1.21.0
  - ggplot2: v3.3.5
  - knitr: v1.36
  - MASS: v7.3-54
  - parallel: v4.1.1
  - rmarkdown: v2.8
  - testthat: v3.0.2
