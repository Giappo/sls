context("sls_dummy")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("logliks", {
  functions <- sls::sls_logliks_dummy()
  testthat::expect_true(length(functions) > 0)
  pars <- c(0.3, 0.2)
  brts <- list(c(10, 9, 6), c(5, 3, 2))
  cond <- 2
  for (i in 1:length(functions)) {
    fun <- eval(functions[i])[[1]]
    test <- fun(
      pars = pars,
      brts = brts,
      cond = cond
    )
    testthat::expect_true(is.numeric(test))
    testthat::expect_true(test <= 0)
  }
})

test_that("ml", {
  if (!is_on_ci()) {
    skip("This has to run on ci.")
  }
  functions <- sls::sls_logliks_dummy()
  testthat::expect_true(length(functions) > 0)
  pars <- c(0.3, 0.2)
  brts <- list(c(10, 9, 6), c(5, 3, 2))
  cond <- 2
  for (i in 1:length(functions)) {
    fun <- eval(functions[i])[[1]]
    test <- sls_ml_dummy(
      loglik_function = fun,
      brts = brts,
      start_pars = c(0.5, 0.3),
      n_0 = 2,
      cond = cond,
      verbose = FALSE
    )
    testthat::expect_true(
      nrow(test) >= 1 && ncol(test) >= 1
    )
  }
})

test_that("main", {
  if (!is_on_ci()) {
    skip("This has to run on ci.")
  }
  seed <- 14
  sim_pars <- c(0.2, 0.1, 0.6, 0.1)
  cond <- 2
  l_2 <- sls::sim_get_standard_l_2(crown_age = 6, shift_time = 4)
  test <- sls::sls_main_dummy(
    seed = seed,
    sim_pars = sim_pars,
    cond = cond,
    l_2 = l_2,
    verbose = FALSE
  )
  testthat::expect_true(
    nrow(test) >= 1 && ncol(test) >= 1
  )
})
