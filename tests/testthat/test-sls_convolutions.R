context("sls_convolutions")

test_that("dft is a unitary operator", {
  vec <- 1:30
  testthat::expect_true(
    all(unlist(Re(sls::idft(sls::dft(vec)))) - vec < 1e-10)
  )
  testthat::expect_true(
    all(unlist(Re(sls::dft(sls::idft(vec)))) - vec < 1e-10)
  )
})

test_that("dft yields the same result as the standard convolution", {
  lambda <- 0.3
  mu <- 0.1
  times <- c(10, 8, 7, 3)
  tbar <- 2
  n_max <- 15
  res_conv <- sls::combine_pns0(
    lambda = lambda,
    mu = mu,
    times = times,
    tbar = tbar,
    n_max = n_max
  )
  res_dft  <- sls::combine_pns(
    lambda = lambda,
    mu = mu,
    times = times,
    tbar = tbar,
    n_max = n_max
  )
  testthat::expect_equal(res_conv, res_dft)
})

test_that("ddd yields the same result as the other methods", {
  for (ii in 1:40) {
    set.seed(ii)
    lambda <- runif(n = 1, min = 0.1, max = 1)
    mu <- runif(n = 1, min = 0, max = lambda * (4 / 5))
    times <- c(10, 8.5, 8, 7, 6, 5.7, 3)
    tbar <- 2
    n_max <- 50
    test_sls <- sls::combine_pns(
      lambda = lambda,
      mu = mu,
      times = times,
      tbar = tbar,
      n_max = n_max
    )
    test_ddd <- sls::combine_ddd(
      lambda = lambda,
      mu = mu,
      times = times,
      tbar = tbar,
      n_max = n_max
    )
    testthat::expect_equal(test_sls, test_ddd)
  }
})

test_that("ddd yields the same wrong result as the other wrong methods", {
  for (ii in 1:40) {
    set.seed(ii)
    lambda <- runif(n = 1, min = 0.1, max = 1)
    mu <- runif(n = 1, min = 0, max = lambda * (4 / 5))
    times <- c(10, 8.5, 8, 7, 6, 5.7, 3)
    tbar <- 2
    n_max <- 50
    test_sls <- sls::combine_pns_nodiv(
      lambda = lambda,
      mu = mu,
      times = times,
      tbar = tbar,
      n_max = n_max
    )
    test_ddd <- sls::combine_ddd_nodiv(
      lambda = lambda,
      mu = mu,
      times = times,
      tbar = tbar,
      n_max = n_max
    )
    testthat::expect_equal(test_sls, test_ddd)
  }
})
