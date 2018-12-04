context("sls_ml")

test_that("use", {

  brts_m <- c(3, 2, 1)
  brts_s <- c(2.5, 1.5)
  brts <- list(brts_m, brts_s)
  start_pars <- c(0.4, 0.05, 0.3, 0.1)
  cond <- 3
  n_0 <- 2

  test <- sls_ml(
    loglik_function = sls::loglik_sls_p,
    brts = brts,
    start_pars = start_pars,
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
    test$df == length(start_pars)
  )
  testthat::expect_true(
    test$conv == 0
  )
})

test_that("abuse", {

  brts_m <- c(3, 2, 1)
  brts_s <- c(2.5, 1.5)
  brts <- list(brts_m, brts_s)
  start_pars <- c(0.4, 0.05, 0.3, 0.1)
  cond <- 3
  n_0 <- 2

  testthat::expect_output(
    sls_ml(
      loglik_function = sls::loglik_sls_p,
      brts = brts,
      start_pars = c(0, 1, 0, 1),
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    ),
    "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values." # nolint
  )
  testthat::expect_error(
    test <- sls_ml(
      loglik_function = sls::loglik_sls_p,
      brts = list(c(), brts[[2]]),
      start_pars = start_pars,
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    ),
    "main clade branching times cannot be an empty vector"
  )
  testthat::expect_error(
    test <- sls_ml(
      loglik_function = sls::loglik_sls_p,
      brts = list(brts[[1]], c()),
      start_pars = start_pars,
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    ),
    "sub clade branching times cannot be an empty vector"
  )
  testthat::expect_error(
    test <- sls_ml(
      loglik_function = sls::loglik_sls_p,
      brts = brts,
      start_pars = start_pars,
      cond = 15,
      n_0 = n_0,
      verbose = FALSE
    ),
    "this conditioning is not implemented"
  )
  testthat::expect_error(
    test <- sls_ml(
      loglik_function = sls::loglik_sls_p,
      brts = brts,
      start_pars = start_pars,
      cond = cond,
      n_0 = 3,
      verbose = FALSE
    ),
    "this n_0 is not implemented"
  )
})
