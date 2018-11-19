context("convolutions")

test_that("DFT is a unitary operator", {
  vec <- 1:30
  testthat::expect_true(
    all(unlist(Re(sls::IDFT(sls::DFT(vec)))) - vec < 1e-10)
  )
  testthat::expect_true(
    all(unlist(Re(sls::DFT(sls::IDFT(vec)))) - vec < 1e-10)
  )
})

test_that("DFT yields the same result as the standard convolution", {
  lambda <- 0.3
  mu <- 0.1
  ts <- c(10, 8, 7, 3)
  tbar <- 2
  nmax <- 15
  res_conv <- sls::combine_pns0(
    lambda = lambda,
    mu = mu,
    ts = ts,
    tbar = tbar,
    nmax = nmax
  )
  res_dft  <- sls::combine_pns(
    lambda = lambda,
    mu = mu,
    ts = ts,
    tbar = tbar,
    nmax = nmax
  )
  testthat::expect_equal(res_conv, res_dft)
})
