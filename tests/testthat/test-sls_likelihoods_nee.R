context("sls_likelihooods_nee")

test_that("1 minus functions work", {
  lambda <- 0.25
  mu <- 0.1
  t <- 6
  expect_equal(
    1 - sls::pt(t = t, lambda = lambda, mu = mu),
    sls::one_minus_pt(t = t, lambda = lambda, mu = mu)
  )
  expect_equal(
    1 - sls::ut(t = t, lambda = lambda, mu = mu),
    sls::one_minus_ut(t = t, lambda = lambda, mu = mu)
  )
  ###
  lambda <- 0.25
  mu <- 0.25
  t <- 7
  expect_equal(
    1 - sls::pt(t = t, lambda = lambda, mu = mu),
    sls::one_minus_pt(t = t, lambda = lambda, mu = mu)
  )
  expect_equal(
    1 - sls::ut(t = t, lambda = lambda, mu = mu),
    sls::one_minus_ut(t = t, lambda = lambda, mu = mu)
  )
})
