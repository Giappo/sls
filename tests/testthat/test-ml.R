context("sls_ml")

test_that("use", {

  brts_m <- c(3, 2, 1)
  brts_s <- c(2.5, 1.5)
  startpars <- c(0.4, 0.05, 0.3, 0.1)
  cond <- 3
  n_0 <- 2

  test <- sls_ml(
    loglik_function = sls::loglik_slsP,
    brts_m = brts_m,
    brts_s = brts_s,
    startpars = startpars,
    cond = cond,
    n_0 = n_0,
    verbose = FALSE
  )

  testthat::expect_true(
    test$lambda_m >= 0
  )
  testthat::expect_true(
    test$mu_m >= 0
  )
  testthat::expect_true(
    test$lambda_s >= 0
  )
  testthat::expect_true(
    test$mu_s >= 0
  )
  testthat::expect_true(
    is.numeric(test$loglik) == TRUE
  )
  testthat::expect_true(
    test$df == length(startpars)
  )
  testthat::expect_true(
    test$conv == 0
  )
})

test_that("abuse", {

  brts_m <- c(3, 2, 1)
  brts_s <- c(2.5, 1.5)
  startpars <- c(0.4, 0.05, 0.3, 0.1)
  cond <- 3
  n_0 <- 2

  testthat::expect_error(
    test <- sls_ml(
      loglik_function = sls::loglik_slsP,
      brts_m = c(),
      brts_s = brts_s,
      startpars = startpars,
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    ),
    "main clade branching times cannot be an empty vector"
  )
  testthat::expect_error(
    test <- sls_ml(
      loglik_function = sls::loglik_slsP,
      brts_m = brts_m,
      brts_s = c(),
      startpars = startpars,
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    ),
    "sub clade branching times cannot be an empty vector"
  )
  testthat::expect_error(
    test <- sls_ml(
      loglik_function = sls::loglik_slsP,
      brts_m = brts_m,
      brts_s = brts_s,
      startpars = c(-1, startpars[2:4]),
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    ),
    "you cannot start from negative parameters"
  )
  testthat::expect_error(
    test <- sls_ml(
      loglik_function = sls::loglik_slsP,
      brts_m = brts_m,
      brts_s = brts_s,
      startpars = startpars,
      cond = 15,
      n_0 = n_0,
      verbose = FALSE
    ),
    "this conditioning is not implemented"
  )
  testthat::expect_error(
    test <- sls_ml(
      loglik_function = sls::loglik_slsP,
      brts_m = brts_m,
      brts_s = brts_s,
      startpars = startpars,
      cond = cond,
      n_0 = 3,
      verbose = FALSE
    ),
    "this n_0 is not implemented"
  )
})
