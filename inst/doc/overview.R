## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = FALSE, echo = FALSE-------------------------------------------
Sys.setenv(LANGUAGE="en")

## ---- echo = FALSE------------------------------------------------------------
library(robust2sls)
library(AER) # for 2sls regression

## -----------------------------------------------------------------------------
param <- generate_param(dx1 = 1, dx2 = 1, dz2 = 1, intercept = TRUE, beta = c(2,-1), sigma = 1, seed = 2)
data_full <- generate_data(parameters = param, n = 1000)$data

# have a look at the data
head(data_full)

## -----------------------------------------------------------------------------
cor(data_full$u, data_full$x2) # correlation for endogenous regressor
cor(data_full$u, data_full$z2) # close to zero correlation for valid instrument
cor(data_full$x2, data_full$z2) # correlation for informativeness

## -----------------------------------------------------------------------------
fml_ols <- y ~ x2
ols <- lm(formula = fml_ols, data = data_full[, 1:3])
print(ols$coefficients)

fml_tsls <- y ~ x2 | z2
tsls <- ivreg(formula = fml_tsls, data = data_full[, c(1:3, 6)])
print(tsls$coefficients)

## -----------------------------------------------------------------------------
# extract the variables we observe
data <- data.frame(data_full[, c(1:3, 6)])

# not iterating the algorithm
robustified_0 <- outlier_detection(data = data, formula = fml_tsls, ref_dist = "normal",
                                   sign_level = 0.01, initial_est = "robustified", iterations = 0)
print(robustified_0)

# iterating algorithm until convergence
robustified_conv <- outlier_detection(data = data, formula = fml_tsls, ref_dist = "normal",
                                      sign_level = 0.01, initial_est = "robustified", 
                                      iterations = "convergence", convergence_criterion = 0)
print(robustified_conv)


## -----------------------------------------------------------------------------
sum(abs(data_full[, "u"]) > 2.58)

## -----------------------------------------------------------------------------
# which observations had large true errors?
large_true <- which(abs(data_full[, "u"]) > 2.58)

# which observations were detected as outliers
# both step 0 and iterated version had same starting point, so same initial classification
large_detected_0 <- which(robustified_conv$type$m0 == 0)
large_detected_3 <- which(robustified_conv$type$m3 == 0)

# how much do the sets of true and detected large errors overlap?
sum(large_detected_0 %in% large_true) # 12 of the 13 detected outliers really had large errors
sum(large_detected_3 %in% large_true) # 13 of the 15 detected outliers really had large errors

# which were wrongly detected as having large errors?
ind <- large_detected_3[which(!(large_detected_3 %in% large_true))]
# what were their true error values?
data_full[ind, "u"] # relatively large values even though technically smaller than cut-off

# which large true errors were missed?
ind2 <- large_true[which(!(large_true %in% large_detected_3))]
# what were their true error values?
data_full[ind2, "u"] # one of them is close to the cut-off

## -----------------------------------------------------------------------------
# not iterating the algorithm
saturated_0 <- outlier_detection(data = data, formula = fml_tsls, ref_dist = "normal",
                                 sign_level = 0.01, initial_est = "saturated", iterations = 0,
                                 split = 0.5)
print(saturated_0)

# iterating algorithm until convergence
saturated_conv <- outlier_detection(data = data, formula = fml_tsls, ref_dist = "normal",
                                    sign_level = 0.01, initial_est = "saturated", split = 0.5,
                                    iterations = "convergence", convergence_criterion = 0)
print(saturated_conv)

# which did it find in final selection?
large_detected_4 <- which(saturated_conv$type$m4 == 0)
identical(large_detected_3, large_detected_4)

## -----------------------------------------------------------------------------
# final model (iteration 3) without corrected standard errors
summary(robustified_conv$model$m3)$coef

# final model with corrected standard errors, m = 3 (subset of output)
beta_inf(robustified_conv, iteration = 3, exact = TRUE)[, 1:5]

# final model with corrected standard errors, m = "fixed point" / "converged" (subset of output)
beta_inf(robustified_conv, iteration = 3, exact = TRUE, fp = TRUE)[, 1:5]

## -----------------------------------------------------------------------------
# use case re-sampling (to save time, use low iteration m = 1)
resampling <- case_resampling(robustified_conv, R = 1000, m = 1)
beta_boot <- evaluate_boot(resampling, iterations = 1)
mat <- matrix(beta_boot[, 1:2], ncol = 1, nrow = 2) # show subset, put into column vector
colnames(mat) <- "bootStd. Error"
cbind(beta_inf(robustified_conv, iteration = 1, exact = TRUE)[, 1:3], mat)

## -----------------------------------------------------------------------------
# t-test on a single coefficient
# testing difference in beta_2 (coef on endogenous regressor x2), m = 3
beta_t(robustified_conv, iteration = 3, element = "x2")
# testing difference in beta_2 (coef on endogenous regressor x2), m = "fixed point"
beta_t(robustified_conv, iteration = 3, element = "x2", fp = TRUE)

# Hausman-type test on whole coefficient vector, m = 3
beta_hausman(robustified_conv, iteration = 3)

## -----------------------------------------------------------------------------
data_full_contaminated <- data_full
outlier_location <- sample(1:NROW(data_full), size = 0.03*NROW(data_full), replace = FALSE)
outlier_size <- sample(c(-3.5, 3.5), size = length(outlier_location), replace = TRUE)

# replace errors
data_full_contaminated[outlier_location, "u"] <- outlier_size
# replace value of dependent variable
data_full_contaminated[outlier_location, "y"] <- data_full_contaminated[outlier_location, "y"] - 
                                                  data_full[outlier_location, "u"] + 
                                                  data_full_contaminated[outlier_location, "u"]

# extract the data that researcher would actually collect
data_cont <- data.frame(data_full_contaminated[, c(1:3, 6)])

## -----------------------------------------------------------------------------
# not iterating the algorithm
robustified_0 <- outlier_detection(data = data_cont, formula = fml_tsls, ref_dist = "normal",
                                   sign_level = 0.01, initial_est = "robustified", iterations = 0)
print(robustified_0)

# iterating algorithm until convergence
robustified_conv <- outlier_detection(data = data_cont, formula = fml_tsls, ref_dist = "normal",
                                      sign_level = 0.01, initial_est = "robustified", 
                                      iterations = "convergence", convergence_criterion = 0)
print(robustified_conv)

## -----------------------------------------------------------------------------
# which observations were classified as outliers?
detected_0 <- which(robustified_conv$type$m0 == 0)
detected_7 <- which(robustified_conv$type$m7 == 0)

# overlap between deterministic outliers
sum(outlier_location %in% detected_0) # initial found all 30 (+1 in addition)
sum(outlier_location %in% detected_7) # converged found all 30 (+ 15 in addition)

## -----------------------------------------------------------------------------
# full sample model
tsls_cont <- ivreg(formula = fml_tsls, data = data_cont)
tsls_cont$coefficients

# updated estimates after initial classification and removal of outliers
robustified_conv$model$m1$coefficients
# updated estimates after convergence
robustified_conv$model$m7$coefficients

## -----------------------------------------------------------------------------
# use case re-sampling (to save time, use low iteration m = 1)
resampling <- case_resampling(robustified_conv, R = 1000, m = 1)
beta_boot <- evaluate_boot(resampling, iterations = 1)
mat <- matrix(beta_boot[, 1:2], ncol = 1, nrow = 2) # show subset, put into column vector
colnames(mat) <- "bootStd. Error"
cbind(beta_inf(robustified_conv, iteration = 1, exact = TRUE)[, 1:3], mat)

## -----------------------------------------------------------------------------
# not iterating the algorithm
saturated_0 <- outlier_detection(data = data_cont, formula = fml_tsls, ref_dist = "normal",
                                 sign_level = 0.01, initial_est = "saturated", iterations = 0,
                                 split = 0.5)
print(saturated_0)

# iterating algorithm until convergence
saturated_conv <- outlier_detection(data = data_cont, formula = fml_tsls, ref_dist = "normal",
                                    sign_level = 0.01, initial_est = "saturated", split = 0.5,
                                    iterations = "convergence", convergence_criterion = 0)
print(saturated_conv)

# which did it find in final selection?
detected_7_sat <- which(saturated_conv$type$m7 == 0)
identical(detected_7, detected_7_sat)

