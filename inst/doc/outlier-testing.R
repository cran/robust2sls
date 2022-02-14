## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = FALSE, echo = FALSE-------------------------------------------
Sys.setenv(LANGUAGE="en")

## ---- data--------------------------------------------------------------------
library(robust2sls)
# create parameters
p <- generate_param(3, 2, 3, sigma = 2, intercept = TRUE, seed = 42)
# draw random sample of 1000 observations following the model
d <- generate_data(parameters = p, n = 1000)$data

## ---- single_alg--------------------------------------------------------------
model <- outlier_detection(data = d, formula = p$setting$formula, ref_dist = "normal",
                           sign_level = 0.05, initial_est = "robustified", iterations = 5)
print(model)

## ---- multiple_alg------------------------------------------------------------
# choose which gamma values to use
gammas <- seq(0.01, 0.05, 0.01)
# register backend
library(doFuture)
registerDoFuture()
plan(sequential)
models <- multi_cutoff(gamma = gammas, data = d, formula = p$setting$formula, ref_dist = "normal",
                       initial_est = "robustified", iterations = 5)
length(models)

## ---- proptest----------------------------------------------------------------
# using a single robust2sls object
proptest(model, alpha = 0.05, iteration = 5, one_sided = FALSE)

# using a list of robust2sls objects
proptest(models, alpha = 0.05, iteration = 5, one_sided = TRUE)

## ---- proptest simes----------------------------------------------------------
proptests <- proptest(models, alpha = 0.05, iteration = 5, one_sided = TRUE)
a <- globaltest(tests = proptests, global_alpha = 0.05)

# decision for global hypothesis test
a$reject

# details for the Simes procedure
a$tests[, c("iter_test", "iter_act", "gamma", "t", "pval", "alpha_adj", "reject_adj")]


## ---- counttest---------------------------------------------------------------
# using a single robust2sls object
counttest(model, alpha = 0.05, iteration = 5, one_sided = FALSE)

# using a list of robust2sls objects
counttest(models, alpha = 0.05, iteration = 5, one_sided = TRUE)

## ---- counttest simes---------------------------------------------------------
counttests <- counttest(models, alpha = 0.05, iteration = 5, one_sided = TRUE)
b <- globaltest(tests = counttests, global_alpha = 0.05)

# decision for global hypothesis test
b$reject

# details for the Simes procedure
b$tests[, c("iter_test", "iter_act", "gamma", "num_act", "num_exp", "pval", "alpha_adj", "reject_adj")]

## ---- sumtest-----------------------------------------------------------------
c <- sumtest(models, alpha = 0.05, iteration = 1, one_sided = FALSE)

attr(c, "gammas")

## ---- suptest-----------------------------------------------------------------
d <- suptest(models, alpha = 0.05, iteration = 5)

attr(d, "gammas")
attr(d, "critical")

