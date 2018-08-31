context("conditional probabilities")

test_that("Analytical equivalence for conditional probability 1", {
  lambdas    <- c(0.5, 0.4)
  mus        <- c(0.4, 0.2)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  tp <- 0 ;tc <- -10; ts <- -6

  ns <- 2:1e5
  test1 <- sum(
    sls::pn(n = ns, t = ts - tc, lambda = lambdas[1], mu = mus[1]) *
      (1 - (sls::pn(n = 0, t = tp - ts, lambda = lambdas[1], mu = mus[1])) ^ (ns - 1) )
  )


  test2 <- sls::pn(n = 1, t = ts - tc, lambda = lambdas[1], mu = mus[1]) *
    (1 - sls::pn(n = 0, t = tp - ts, lambda = lambdas[1], mu = mus[1])) *
    sls::ut(t = ts - tc, lambda = lambdas[1], mu = mus[1]) *
    (1 - sls::ut(t = ts - tc, lambda = lambdas[1], mu = mus[1]))^-1 *
    (1 - sls::ut(t = ts - tc, lambda = lambdas[1], mu = mus[1]) * sls:::pn(n = 0, t = tp - ts, lambda = lambdas[1], mu = mus[1]))^-1

  testthat::expect_equal(
    test1, test2
  )

})

test_that("Pc1 never smaller than Pc3 (Pc3 is a stricter condition)", {

  #test1
  brtsM      <- c(-12, -9)
  brtsS      <- c(-5)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  conditioning <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_true(conditioning[[1]] >= conditioning[[3]])

  #test2
  brtsM      <- c(-12, -9, -4)
  brtsS      <- c(-5)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  conditioning <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_true(conditioning[[1]] >= conditioning[[3]])

  #test3
  brtsM      <- c(-12, -9, -4)
  brtsS      <- c(-5)
  lambdas    <- c(0.5, 0.7)
  mus        <- c(0.4, 0.2)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  conditioning <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_true(conditioning[[1]] >= conditioning[[3]])

  #test4
  brtsM      <- c(-12, -9, -4, -2)
  brtsS      <- c(-5)
  lambdas    <- c(0.5, 0.4)
  mus        <- c(0.4, 0.3)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  conditioning <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_true(conditioning[[1]] >= conditioning[[3]])

  #test5
  brtsM      <- c(-10, -10, -9, -6, -4, -1)
  brtsS      <- c(-5, -3)
  lambdas    <- c(0.4, 0.6)
  mus        <- c(0.3, 0.1)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  conditioning <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_true(conditioning[[1]] >= conditioning[[3]])
})

test_that("If lambda2=mu2=0 (inert subclade), Pc1 equal to Pc3 (PS = 1)", {

  brtsM      <- c(-10, -10, -9, -6, -4, -1)
  brtsS      <- c(-5, -3)
  lambdas    <- c(0.4, 0)
  mus        <- c(0.3, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  conditioning <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_true(conditioning[[1]] == conditioning[[3]])
})
