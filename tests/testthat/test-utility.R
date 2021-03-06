test_that("extract_formula() works correctly", {

  f1 <- y ~ x1 + x2 | x1 + z2 + z3
  v1 <- extract_formula(f1)

  expect_equal(length(v1),5)
  expect_equal(length(v1$y_var), 1)
  expect_equal(length(v1$x1_var), 1)
  expect_equal(length(v1$x2_var), 1)
  expect_equal(length(v1$z1_var), 1)
  expect_equal(length(v1$z2_var), 2)
  expect_equal(v1$y_var, "y")
  expect_equal(v1$x1_var, "x1")
  expect_equal(v1$x2_var, "x2")
  expect_equal(v1$z1_var, "x1")
  expect_equal(v1$z2_var, c("z2", "z3"))

  expect_error(extract_formula(y ~ x1 + x2 + z2), "does not include both symbols
    `~` and `|`")
  expect_error(extract_formula(y ~ z1 + z2), "does not include both symbols
    `~` and `|`")
  expect_error(extract_formula(~ x1 + x2 | x1 + z2), "any dependent variable")
  # expect_error(extract_formula(y ~ | z1 + z2)) # this is not possible anyways
  # expect_error(extract_formula(y ~ x1 + x2 |)) # this is not possible anyways

  f2 <- mpg ~ cyl * disp | cyl + wt + gear
  v2 <- extract_formula(f2)

  expect_equal(v2$y_var, "mpg")
  expect_equal(v2$x1_var, "cyl")
  expect_equal(v2$x2_var, "disp")
  expect_equal(v2$z1_var, "cyl")
  expect_equal(v2$z2_var, c("wt", "gear"))

  f3 <- y ~ x1 + x2 + x3 | x1 + z2
  expect_error(extract_formula(f3), "does not fulfill the order condition")

})

test_that("selection() works correctly", {
  dta <- datasets::mtcars
  dta[1, "mpg"] <- NA
  dta[2, "cyl"] <- NA
  dta[3, "wt"] <- NA
  test <- AER::ivreg(mpg ~ cyl + disp | cyl + wt, data = dta)

  expect_equal(NROW(dta), 32) # 32 observations
  expect_equal(NROW(test$residuals), 29) # 29 obs due to NA values in y, x, z

  sel <- selection(data = dta, yvar = "mpg", model = test, cutoff = 1.96)
  res1 <- res2 <- res3 <- as.double(NA)
  names(res1) <- "Mazda RX4"
  names(res2) <- "Mazda RX4 Wag"
  names(res3) <- "Datsun 710"

  expect_equal(NROW(sel[["res"]]), 32) # 32 when done manually inside res_all()
  # check that these residuals with missing y, x, z are NA
  expect_equal(sel[["res"]][1], res1)
  expect_equal(sel[["res"]][2], res2)
  expect_equal(sel[["res"]][3], res3)
  # check that the model residuals coincide with the res_all(), ignoring NAs
  expect_equal(sel[["res"]][c(-1,-2,-3)], test$residuals)

  # standardised residuals of first three should also be NA
  expect_equal(sel[["stdres"]][1], res1)
  expect_equal(sel[["stdres"]][2], res2)
  expect_equal(sel[["stdres"]][3], res3)

  # in selection vector obs with missing values should be set to FALSE
  sel1 <- sel2 <- sel3 <- FALSE
  names(sel1) <- "Mazda RX4"
  names(sel2) <- "Mazda RX4 Wag"
  names(sel3) <- "Datsun 710"
  expect_equal(sel[["sel"]][1], sel1)
  expect_equal(sel[["sel"]][2], sel2)
  expect_equal(sel[["sel"]][3], sel3)

  # by inspection, only one observation has an absolute residual > 1.96
  # this is "Pontiac Firebird", observation 25 in the original dataframe
  sel25 <- FALSE
  names(sel25) <- "Pontiac Firebird"
  expect_equal(sum(sel[["sel"]]), 28) # 28 because 3 missing and one outlier
  expect_equal(sel[["sel"]][c(1, 2, 3, 25)], c(sel1, sel2, sel3, sel25))
  expect_equal(as.vector(sel[["sel"]][c(-1, -2, -3, -25)]),
               !logical(length = 28)) # strip the names attribute via as.vector

  # check that type has the correct coding
  expect_equal(typeof(sel[["type"]]), "integer")
  expect_equal(as.vector(sel[["type"]][1]), -1)
  expect_equal(as.vector(sel[["type"]][2]), -1)
  expect_equal(as.vector(sel[["type"]][3]), -1)
  expect_equal(as.vector(sel[["type"]][4]), 1)
  expect_equal(as.vector(sel[["type"]][25]), 0)
  expect_equal(as.vector(sel[["type"]][c(1,2,3)]), c(-1, -1, -1))

  test <- AER::ivreg(mpg ~ cyl + disp | cyl + wt, data = mtcars)
  sel <- selection(data = mtcars, yvar = "mpg", model = test, cutoff = 1.96)
  expect_equal(sel[["res"]], test$residuals)

  lmmodel <- stats::lm(mpg ~ cyl + disp, data = mtcars)
  expect_error(selection(data = mtcars, yvar = "mpg", model = lmmodel,
                         cutoff = 1.96), "not of class `ivreg`")
  expect_error(selection(data = mtcars, yvar = 1, model = test, cutoff = 1.96),
               "not a character vector")
  expect_error(selection(data = mtcars, yvar = "notexist", model = test,
                         cutoff = 1.96), "not a variable in the dataframe")




})

test_that("nonmissing() works correctly", {

  expect_equal(nonmissing(data = mtcars, formula = mpg ~ cyl + disp | cyl + wt),
               !logical(length = 32))

  dta <- datasets::mtcars
  dta[1, "mpg"] <- NA
  dta[2, "cyl"] <- NA
  dta[3, "wt"] <- NA
  dta[4:5, "disp"] <- NA
  dta[15, "gear"] <- NA # gear NA in obs 15 but not used in estimation
  dta[16, "qsec"] <- NA # qsec NA in obs 16 but not used in estimation

  nm <- !logical(length = 32)
  nm[1:5] <- FALSE
  expect_equal(nonmissing(data = dta, formula = mpg ~ cyl + disp | cyl + wt),
               nm)

  expect_error(nonmissing(data = as.matrix(dta),
      formula = mpg ~ cyl + disp | cyl + wt), "must be a dataframe")
  expect_error(nonmissing(data = dta, formula = mpg ~ cyl + disp), "does not
               include both symbols `~` and `|`")
  expect_error(nonmissing(data = dta, formula = ~ cyl | disp), "does not
               include both symbols `~` and `|`")
  expect_error(nonmissing(data = dta, formula = nonexist ~ cyl + disp))

})

test_that("constants() works correctly", {

  # working values
  data <- datasets::mtcars
  # since formula has an environment, whose memory address changes each time
  # it is run, it differs by snapshot. So here remove environment.
  formula <- mpg ~ cyl + disp | cyl + wt
  attr(formula, which = ".Environment") <- NULL
  # attention: captures wrong call (namely of testing environment)
  test1 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "robustified",
            iterations = 5, convergence_criterion = NULL, shuffle = FALSE,
            shuffle_seed = 42, split = 0.5, max_iter = NULL)
  call1 <- sys.call()
  test2 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "robustified",
            iterations = 5, convergence_criterion = 0, shuffle = NULL,
            shuffle_seed = NULL, split = NULL, max_iter = NULL)
  call2 <- sys.call()
  test3 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "saturated",
            iterations = "convergence", convergence_criterion = 0.5,
            shuffle = TRUE, shuffle_seed = 42, split = 0.5, max_iter = NULL)
  call3 <- sys.call()
  test4 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "saturated",
            iterations = "convergence", convergence_criterion = 1,
            shuffle = FALSE, shuffle_seed = 42, split = 0.5, max_iter = NULL)
  call4 <- sys.call()

  c1 <- constants(call = call1, formula = formula, data = data,
                  reference = "normal", sign_level = 0.05,
                  estimator = "robustified", split = 0.5, shuffle = FALSE,
                  shuffle_seed = 42, iter = 5, criterion = NULL,
                  user_model = NULL, verbose = FALSE, max_iter = NULL)
  c2 <- constants(call = call2, formula = formula, data = data,
                  reference = "normal", sign_level = 0.05,
                  estimator = "robustified", split = NULL, shuffle = NULL,
                  shuffle_seed = NULL, iter = 5, criterion = 0,
                  user_model = NULL, verbose = FALSE, max_iter = NULL)
  c3 <- constants(call = call3, formula = formula, data = data,
                  reference = "normal", sign_level = 0.05,
                  estimator = "saturated", split = 0.5, shuffle = TRUE,
                  shuffle_seed = 42, iter = 5, criterion = 0.5,
                  user_model = NULL, verbose = FALSE, max_iter = NULL)
  c4 <- constants(call = call4, formula = formula, data = data,
                  reference = "normal", sign_level = 0.05,
                  estimator = "saturated", split = 0.5, shuffle = FALSE,
                  shuffle_seed = 42, iter = "convergence", criterion = 1,
                  user_model = NULL, verbose = FALSE, max_iter = NULL)
  c5 <- constants(call = call4, formula = formula, data = data,
                  reference = "normal", sign_level = 0.05,
                  estimator = "saturated", split = 0.5, shuffle = FALSE,
                  shuffle_seed = 42, iter = "convergence", criterion = 1,
                  user_model = NULL, verbose = FALSE, max_iter = 10)

  names <- c("call", "verbose", "formula", "data", "reference", "sign_level",
             "psi", "cutoff", "bias_corr", "initial", "convergence",
             "iterations")

  expect_equal(length(c1), 12)
  expect_equal(length(c2), 12)
  expect_equal(length(c3), 12)
  expect_equal(length(c4), 12)
  expect_equal(length(c5), 12)
  expect_equal(class(c1), "list")
  expect_equal(class(c2), "list")
  expect_equal(class(c3), "list")
  expect_equal(class(c4), "list")
  expect_equal(class(c5), "list")
  expect_equal(names(c1), names)
  expect_equal(names(c2), names)
  expect_equal(names(c3), names)
  expect_equal(names(c4), names)
  expect_equal(names(c5), names)
  expect_equal(length(c1$initial), 5)
  expect_equal(length(c2$initial), 5)
  expect_equal(length(c3$initial), 5)
  expect_equal(length(c4$initial), 5)
  expect_equal(length(c5$initial), 5)
  expect_equal(class(c1$initial), "list")
  expect_equal(class(c2$initial), "list")
  expect_equal(class(c3$initial), "list")
  expect_equal(class(c4$initial), "list")
  expect_equal(class(c5$initial), "list")
  expect_equal(names(c1$initial),
               c("estimator", "split", "shuffle", "shuffle_seed", "user"))
  expect_equal(names(c2$initial),
               c("estimator", "split", "shuffle", "shuffle_seed", "user"))
  expect_equal(names(c3$initial),
               c("estimator", "split", "shuffle", "shuffle_seed", "user"))
  expect_equal(names(c4$initial),
               c("estimator", "split", "shuffle", "shuffle_seed", "user"))
  expect_equal(names(c5$initial),
               c("estimator", "split", "shuffle", "shuffle_seed", "user"))
  expect_equal(length(c1$convergence), 5)
  expect_equal(length(c2$convergence), 5)
  expect_equal(length(c3$convergence), 5)
  expect_equal(length(c4$convergence), 5)
  expect_equal(length(c5$convergence), 5)
  expect_equal(class(c1$convergence), "list")
  expect_equal(class(c2$convergence), "list")
  expect_equal(class(c3$convergence), "list")
  expect_equal(class(c4$convergence), "list")
  expect_equal(class(c5$convergence), "list")
  expect_equal(names(c1$convergence), c("criterion", "difference", "converged", "iter", "max_iter"))
  expect_equal(names(c2$convergence), c("criterion", "difference", "converged", "iter", "max_iter"))
  expect_equal(names(c3$convergence), c("criterion", "difference", "converged", "iter", "max_iter"))
  expect_equal(names(c4$convergence), c("criterion", "difference", "converged", "iter", "max_iter"))
  expect_equal(names(c5$convergence), c("criterion", "difference", "converged", "iter", "max_iter"))
  expect_equal(length(c1$iterations), 2)
  expect_equal(length(c2$iterations), 2)
  expect_equal(length(c3$iterations), 2)
  expect_equal(length(c4$iterations), 2)
  expect_equal(length(c5$iterations), 2)
  expect_equal(class(c1$iterations), "list")
  expect_equal(class(c2$iterations), "list")
  expect_equal(class(c3$iterations), "list")
  expect_equal(class(c4$iterations), "list")
  expect_equal(class(c5$iterations), "list")
  expect_equal(names(c1$iterations), c("setting", "actual"))
  expect_equal(names(c2$iterations), c("setting", "actual"))
  expect_equal(names(c3$iterations), c("setting", "actual"))
  expect_equal(names(c4$iterations), c("setting", "actual"))
  expect_equal(names(c5$iterations), c("setting", "actual"))

  expect_snapshot_output(c1)
  expect_snapshot_output(c2)
  expect_snapshot_output(c3)
  expect_snapshot_output(c4)
  expect_snapshot_output(c5)


  expect_equal(c1$formula, formula)
  expect_equal(c2$formula, formula)
  expect_equal(c3$formula, formula)
  expect_equal(c4$formula, formula)
  expect_equal(c5$formula, formula)

  expect_equal(c1$verbose, FALSE)
  expect_equal(c2$verbose, FALSE)
  expect_equal(c3$verbose, FALSE)
  expect_equal(c4$verbose, FALSE)
  expect_equal(c5$verbose, FALSE)

  expect_equal(c1$data, mtcars)
  expect_equal(c2$data, mtcars)
  expect_equal(c3$data, mtcars)
  expect_equal(c4$data, mtcars)
  expect_equal(c5$data, mtcars)

  expect_equal(c1$reference, "normal")
  expect_equal(c2$reference, "normal")
  expect_equal(c3$reference, "normal")
  expect_equal(c4$reference, "normal")
  expect_equal(c5$reference, "normal")

  expect_equal(c1$sign_level, 0.05)
  expect_equal(c2$sign_level, 0.05)
  expect_equal(c3$sign_level, 0.05)
  expect_equal(c4$sign_level, 0.05)
  expect_equal(c5$sign_level, 0.05)

  expect_equal(c1$psi, 0.95)
  expect_equal(c2$psi, 0.95)
  expect_equal(c3$psi, 0.95)
  expect_equal(c4$psi, 0.95)
  expect_equal(c5$psi, 0.95)

  expect_equal(c1$cutoff, 1.959964)
  expect_equal(c2$cutoff, 1.959964)
  expect_equal(c3$cutoff, 1.959964)
  expect_equal(c4$cutoff, 1.959964)
  expect_equal(c5$cutoff, 1.959964)

  expect_equal(c1$bias_corr, 1.317798, tolerance = 0.0000001)
  expect_equal(c2$bias_corr, 1.317798, tolerance = 0.0000001)
  expect_equal(c3$bias_corr, 1.317798, tolerance = 0.0000001)
  expect_equal(c4$bias_corr, 1.317798, tolerance = 0.0000001)
  expect_equal(c5$bias_corr, 1.317798, tolerance = 0.0000001)

  expect_equal(c1$initial$estimator, "robustified")
  expect_equal(c2$initial$estimator, "robustified")
  expect_equal(c3$initial$estimator, "saturated")
  expect_equal(c4$initial$estimator, "saturated")
  expect_equal(c5$initial$estimator, "saturated")

  expect_equal(c1$initial$split, NULL)
  expect_equal(c2$initial$split, NULL)
  expect_equal(c3$initial$split, 0.5)
  expect_equal(c4$initial$split, 0.5)
  expect_equal(c5$initial$split, 0.5)

  expect_equal(c1$initial$shuffle, NULL)
  expect_equal(c2$initial$shuffle, NULL)
  expect_equal(c3$initial$shuffle, TRUE)
  expect_equal(c4$initial$shuffle, FALSE)
  expect_equal(c5$initial$shuffle, FALSE)

  expect_equal(c1$initial$shuffle_seed, NULL)
  expect_equal(c2$initial$shuffle_seed, NULL)
  expect_equal(c3$initial$shuffle_seed, 42)
  expect_equal(c4$initial$shuffle_seed, NULL)
  expect_equal(c5$initial$shuffle_seed, NULL)

  expect_equal(c1$initial$user, NULL)
  expect_equal(c2$initial$user, NULL)
  expect_equal(c3$initial$user, NULL)
  expect_equal(c4$initial$user, NULL)
  expect_equal(c5$initial$user, NULL)

  expect_equal(c1$iterations$setting, 5)
  expect_equal(c2$iterations$setting, 5)
  expect_equal(c3$iterations$setting, 5)
  expect_equal(c4$iterations$setting, "convergence")
  expect_equal(c5$iterations$setting, "convergence")

  expect_equal(c1$iterations$actual, NULL)
  expect_equal(c2$iterations$actual, NULL)
  expect_equal(c3$iterations$actual, NULL)
  expect_equal(c4$iterations$actual, NULL)
  expect_equal(c5$iterations$actual, NULL)

  expect_equal(c1$convergence$difference, NULL)
  expect_equal(c2$convergence$difference, NULL)
  expect_equal(c3$convergence$difference, NULL)
  expect_equal(c4$convergence$difference, NULL)
  expect_equal(c5$convergence$difference, NULL)

  expect_equal(c1$convergence$converged, NULL)
  expect_equal(c2$convergence$converged, NULL)
  expect_equal(c3$convergence$converged, NULL)
  expect_equal(c4$convergence$converged, NULL)
  expect_equal(c5$convergence$converged, NULL)

  expect_equal(c1$convergence$iter, NULL)
  expect_equal(c2$convergence$iter, NULL)
  expect_equal(c3$convergence$iter, NULL)
  expect_equal(c4$convergence$iter, NULL)
  expect_equal(c5$convergence$iter, NULL)

  expect_equal(c1$convergence$max_iter, NULL)
  expect_equal(c2$convergence$max_iter, NULL)
  expect_equal(c3$convergence$max_iter, NULL)
  expect_equal(c4$convergence$max_iter, NULL)
  expect_equal(c5$convergence$max_iter, 10)

  # test error messages

  expect_error(constants(call = call1, formula = formula, data = data,
      reference = "nonexist", sign_level = 0.05,
      estimator = "robustified", split = 0.5, shuffle = FALSE,
      shuffle_seed = 42, iter = 5, criterion = NULL))
  expect_error(constants(call = call1, formula = formula, data = data,
      reference = "normal", sign_level = 1.1,
      estimator = "robustified", split = 0.5, shuffle = FALSE,
      shuffle_seed = 42, iter = 5, criterion = NULL), "has to be > 0 and < 1")
  expect_error(constants(call = call1, formula = formula, data = data,
      reference = "normal", sign_level = -0.05,
      estimator = "robustified", split = 0.5, shuffle = FALSE,
      shuffle_seed = 42, iter = 5, criterion = NULL), "has to be > 0 and < 1")

})

test_that("update_list() works correctly", {

  # these tests will only test the mechanics and not its actual usage
  # but I emulate the mechanics:
  # start with a list that has the components that will be updated + another one
  # then update using a list that has the components
  # check whether updated version is correct

  existing <- list(other = 1:10, model = NULL, res = NULL, stdres = NULL,
                   sel = NULL, type = NULL)
  new <- list(model = "a", res = "b", stdres = "c", sel = "d", type = "e")

  u0 <- update_list(current_list = existing, new_info = new, name = "m0")
  u1 <- update_list(current_list = existing, new_info = new, name = "m1")
  u2 <- update_list(current_list = u0, new_info = new, name = "m1")
  u2 <- update_list(current_list = u2, new_info = new, name = "m2")

  expect_equal(u0$other, 1:10)
  expect_equal(u0$model$m0, "a")
  expect_equal(u0$res$m0, "b")
  expect_equal(u0$stdres$m0, "c")
  expect_equal(u0$sel$m0, "d")
  expect_equal(u0$type$m0, "e")

  expect_equal(u1$other, 1:10)
  expect_equal(u1$model$m1, "a")
  expect_equal(u1$res$m1, "b")
  expect_equal(u1$stdres$m1, "c")
  expect_equal(u1$sel$m1, "d")
  expect_equal(u1$type$m1, "e")

  expect_equal(length(u2$other), 10)
  expect_equal(length(u2$model), 3)
  expect_equal(length(u2$res), 3)
  expect_equal(length(u2$stdres), 3)
  expect_equal(length(u2$sel), 3)
  expect_equal(length(u2$type), 3)
  expect_equal(u2$other, 1:10)
  expect_equal(u2$model$m0, "a")
  expect_equal(u2$res$m0, "b")
  expect_equal(u2$stdres$m0, "c")
  expect_equal(u2$sel$m0, "d")
  expect_equal(u2$type$m0, "e")
  expect_equal(u2$model$m1, "a")
  expect_equal(u2$res$m1, "b")
  expect_equal(u2$stdres$m1, "c")
  expect_equal(u2$sel$m1, "d")
  expect_equal(u2$type$m1, "e")
  expect_equal(u2$model$m2, "a")
  expect_equal(u2$res$m2, "b")
  expect_equal(u2$stdres$m2, "c")
  expect_equal(u2$sel$m2, "d")
  expect_equal(u2$type$m2, "e")

  expect_snapshot_output(u0)
  expect_snapshot_output(u1)
  expect_snapshot_output(u2)

})

test_that("conv_difference() works correctly", {

  # working values
  data <- datasets::mtcars
  # since formula has an environment, whose memory address changes each time
  # it is run, it differs by snapshot. So here remove environment.
  formula <- mpg ~ cyl + disp | cyl + wt
  attr(formula, which = ".Environment") <- NULL

  test1 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "robustified",
            iterations = 5, convergence_criterion = NULL, shuffle = FALSE,
            shuffle_seed = 42, split = 0.5)
  test2 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "robustified",
            iterations = 5, convergence_criterion = 0, shuffle = NULL,
            shuffle_seed = NULL, split = NULL)
  test3 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "saturated",
            iterations = "convergence", convergence_criterion = 0.5,
            shuffle = TRUE, shuffle_seed = 42, split = 0.5)
  test4 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "saturated",
            iterations = "convergence", convergence_criterion = 1,
            shuffle = FALSE, shuffle_seed = 42, split = 0.5)

  # notice that irrespective of counter >= 1, always calculate most recent diff
  # I calculated the values manually from the objects and inserted these as the
  # 'expected' values for the tests

  expect_equal(conv_diff(current = test1, counter = 1), 0)
  expect_equal(conv_diff(current = test2, counter = 1), 0)
  expect_equal(conv_diff(current = test3, counter = 3), 0)
  # when counter == 1 and initial == 'saturated' then calculate L2 norm between
  # m1 and both initial m0 models -> get larger one
  # sum((test3$model$m1$coefficients-test3$model$m0$split2$coefficients)^2)
  # sum((test3$model$m1$coefficients-test3$model$m0$split1$coefficients)^2)
  expect_equal(conv_diff(current = test3, counter = 1), 7.269473,
               tolerance = 0.000001)
  expect_equal(conv_diff(current = test4, counter = 3), 0)
  expect_equal(conv_diff(current = test4, counter = 1), 14.18946,
               tolerance = 0.000001)

  test5 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "saturated",
            iterations = 0, convergence_criterion = NULL,
            shuffle = FALSE, shuffle_seed = 42, split = 0.5)
  test6 <- outlier_detection(data = data, formula = formula,
            ref_dist = "normal", sign_level = 0.05, initial_est = "robustified",
            iterations = 0, convergence_criterion = 0,
            shuffle = FALSE, shuffle_seed = 42, split = 0.5)
  # when only estimate initial model then cannot calculate a difference
  # function should return NULL (value of counter doesn't matter in this case)
  expect_equal(conv_diff(current = test5, counter = 1), NULL)
  expect_equal(conv_diff(current = test5, counter = 3), NULL)
  expect_equal(conv_diff(current = test6, counter = 1), NULL)
  expect_equal(conv_diff(current = test6, counter = 3), NULL)

})



test_that("varrho() works correctly", {

  # set some parameters
  sign_level <- 0.05
  ref_dist = "normal"
  iteration = 2

  # calculate the terms manually then compare to function result
  gamma <- sign_level
  c <- stats::qnorm(gamma/2, mean=0, sd=1, lower.tail = FALSE)
  phi <- 1 - gamma
  f <- stats::dnorm(c, mean=0, sd=1)
  tau_c_2 <- phi - 2 * c * f
  tau_c_4 <- 3 * phi - 2 * c * (c^2 + 3) * f
  tau_2 <- 1
  tau_4 <- 3
  varsigma_c_2 <- tau_c_2 / phi

  vbb <- (2 * c * f / phi)^iteration
  vss <- (c * (c^2 - varsigma_c_2) * f / tau_c_2)^iteration
  vbxu <- (phi^iteration - (2 * c * f)^iteration) /
    (phi^iteration * (phi - 2 * c * f))
  vsuu <- (tau_c_2^iteration - (c*(c^2-varsigma_c_2)*f)^iteration) /
    (tau_c_2^iteration*(tau_c_2 - c*(c^2 - varsigma_c_2)*f))

  vsb0 <- (2*c*f/phi) * 1
  vsb1 <- 1 * (c * (c^2 - varsigma_c_2) * f / tau_c_2)
  vsb <- vsb0 + vsb1

  vsxu0 <- (2*c*f/phi) * 1
  vsxu1 <- 1 * c * (c^2 - varsigma_c_2) * f / tau_c_2
  vsxuterm1 <- (tau_c_2^iteration - (c*(c^2 - varsigma_c_2)*f)^iteration) /
    (tau_c_2^(iteration-1) * (tau_c_2 - c*(c^2 - varsigma_c_2)*f))
  vsxuterm2 <- vsxu0 + vsxu1
  vsxuterm3 <- f / (tau_c_2*(phi - 2*c*f))
  vsxu <- (vsxuterm1 - vsxuterm2) * vsxuterm3

  res <- varrho(sign_level = sign_level, ref_dist = ref_dist,
                iteration = iteration)

  expect_equal(res$c$vbb, vbb)
  expect_equal(res$c$vss, vss)
  expect_equal(res$c$vbxu, vbxu)
  expect_equal(res$c$vsuu, vsuu)
  expect_equal(res$c$vsb, vsb)
  expect_equal(res$c$vsxu, vsxu)

  # test correct error messages
  expect_error(varrho("a", "normal", 1),
               "'sign_level' must be a numeric vector of length 1")
  expect_error(varrho(c(0.01, 0.05), "normal", 1),
               "'sign_level' must be a numeric vector of length 1")
  expect_error(varrho(1.3, "normal", 1),
               "'sign_level' must lie strictly between 0 and 1")
  expect_error(varrho(-0.05, "normal", 1),
               "'sign_level' must lie strictly between 0 and 1")
  expect_error(varrho(0.01, 1, 1),
               "'ref_dist' must be a character vector of length 1")
  expect_error(varrho(0.01, c("normal", "normal"), 1),
               "'ref_dist' must be a character vector of length 1")
  expect_error(varrho(0.05, "nonexistent", 1),
               "'ref_dist' must be one of the available reference")
  expect_error(varrho(0.05, "normal", "3"),
               "'iteration' must be a numeric vector of length 1")
  expect_error(varrho(0.05, "normal", c(1,5)),
               "'iteration' must be a numeric vector of length 1")
  expect_error(varrho(0.05, "normal", 1.2),
               "'iteration' must be an integer")
  expect_error(varrho(0.05, "normal", 0),
               "'iteration' must be weakly larger than 1")
  expect_error(varrho(0.05, "normal", -1),
               "'iteration' must be weakly larger than 1")

})

test_that("mvn_sup() works correctly", {

  m <- c(0, 1)
  S <- matrix(c(2, 1, 1, 3), 2, 2)

  expect_error(mvn_sup(1.3, m, S), "Argument 'n' is not an integer")
  expect_error(mvn_sup("1", m, S), "Argument 'n' is not an integer")
  expect_error(mvn_sup(10, c("a", "b"), S), "'mu' is not a numeric vector")
  expect_error(mvn_sup(10, m, as.character(S)), "'Sigma' is not a numeric matrix")
  expect_error(mvn_sup(10, m, as.data.frame(S)), "'Sigma' is not a numeric matrix")
  expect_error(mvn_sup(10, c(0, 1, 0), S), "'mu' is not compatible with var-cov matrix 'Sigma'")

  a <- round(mvn_sup(2, m, S, seed = 40), 2)
  b <- mvn_sup(10, m, S, seed = 40)
  expect_identical(a, c(1.34, 1.33))
  expect_snapshot_output(b)

  # all output must be weakly positive because we take absolute value
  c <- mvn_sup(10000, m, S, seed = 40)
  expect_equal((c >= 0), rep(TRUE, 10000))
  # output must have length of n (not mu)
  expect_equal(length(c), 10000L)

})
