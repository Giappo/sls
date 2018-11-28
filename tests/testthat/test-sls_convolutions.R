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
