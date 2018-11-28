context("sls_utils")

test_that("p_transition_matrix", {
  lambda <- 0.5
  mu <- 0.1
  matrix_size <- 20
  test <- p_transition_matrix(
    lambda = lambda,
    mu = mu,
    matrix_size = matrix_size
  )
  testthat::expect_true(
    all(dim(test) == rep(matrix_size, 2))
  )
})

test_that("cat2", {
  testthat::expect_output(
    cat2(
      message = "ciao",
      verbose = TRUE
    )
  )
  testthat::expect_silent(
    cat2(
      message = "ciao",
      verbose = FALSE
    )
  )
})

test_that("sls_check_input", {
  brts_m <- c(6, 3, 2)
  brts_s <- c(2.5, 1)
  start_pars <- c(0.5, 0.3, 0.5, 0.3)
  cond <- 3
  n_0 <- 2
  n_max <- 1e2
  testthat::expect_silent(
    sls_check_input(
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_0 = n_0,
      n_max = n_max
    )
  )
  testthat::expect_error(
    test <- sls_check_input(
      brts_m = c(),
      brts_s = brts_s,
      cond = cond,
      n_0 = n_0,
      n_max = n_max
    ),
    "main clade branching times cannot be an empty vector"
  )
  testthat::expect_error(
    test <- sls_check_input(
      brts_m = brts_m,
      brts_s = c(),
      cond = cond,
      n_0 = n_0,
      n_max = n_max
    ),
    "sub clade branching times cannot be an empty vector"
  )
  testthat::expect_error(
    test <- sls_check_input(
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_0 = n_0,
      n_max = -1
    ),
    "it's not going to work with maximum species set to 0 or less"
  )
  testthat::expect_error(
    test <- sls_check_input(
      brts_m = brts_m,
      brts_s = brts_s,
      cond = 17,
      n_0 = n_0,
      n_max = n_max
    ),
    "this conditioning is not implemented"
  )
  testthat::expect_error(
    test <- sls_check_input(
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_0 = 10,
      n_max = n_max
    ),
    "this n_0 is not implemented"
  )
})

test_that("sls_conds", {
  testthat::expect_true(
    is.numeric(sls_conds())
  )
})

test_that("sls_n_0s", {
  testthat::expect_true(
    is.numeric(sls_n_0s())
  )
})

test_that("sls_logliks_div", {
  testthat::expect_true(
    length(sls_logliks_div()) > 0
  )
})

test_that("sls_logliks_nodiv", {
  testthat::expect_true(
    length(sls_logliks_nodiv()) > 0
  )
})

test_that("sls_pkg_name", {
  testthat::expect_true(
    sls_pkg_name() == "sls"
  )
})
