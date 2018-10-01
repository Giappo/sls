context("likelihoods")

test_that( "test P and Q approach equivalence", {

  #test1
  lambdas <- c(0.3, 0)
  mus     <- c(0.1, 0)
  brtsM   <- c(10, 8, 7, 2)
  brtsS   <- NULL
  tsplit  <- c(7)
  td      <- c(5.5)
  pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
  pars2   <- c(100, 1, 1, tsplit, 0, 2)

  lik_P <- sls::loglik_slsP(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_P
  lik_Q <- sls::loglik_slsQ(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_Q

  testthat::expect_equal(
    lik_P, lik_Q
  )

  #test2
  lambdas <- c(0.3, 0.6)
  mus     <- c(0.1, 0.05)
  brtsM   <- c(10, 8, 7, 2)
  brtsS   <- NULL
  tsplit  <- c(7)
  td      <- c(5.5)
  pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
  pars2   <- c(100, 1, 1, tsplit, 0, 2)

  lik_P <- sls::loglik_slsP(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_P
  lik_Q <- sls::loglik_slsQ(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_Q

  testthat::expect_equal(
    lik_P, lik_Q
  )

  #test3
  lambdas <- c(0.3, 0)
  mus     <- c(0.1, 0)
  brtsM   <- c(10, 8, 7, 2)
  brtsS   <- c(5, 3)
  tsplit  <- c(7)
  td      <- c(5.5)
  pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
  pars2   <- c(100, 1, 1, tsplit, 0, 2)

  lik_P <- sls::loglik_slsP(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_P
  lik_Q <- sls::loglik_slsQ(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_Q

  testthat::expect_equal(
    lik_P, lik_Q
  )

  #test4
  lambdas <- c(0.3, 0.6)
  mus     <- c(0.1, 0.05)
  brtsM   <- c(10, 8, 7, 2)
  brtsS   <- c(5, 3)
  tsplit  <- c(7)
  td      <- c(5.5)
  pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
  pars2   <- c(100, 1, 1, tsplit, 0, 2)

  lik_P <- sls::loglik_slsP(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_P
  lik_Q <- sls::loglik_slsQ(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_Q

  testthat::expect_equal(
    lik_P, lik_Q
  )

  #test5
  lambdas <- c(0.2, 0)
  mus     <- c(0.1, 0)
  brtsM   <- c(10, 9.5, 8, 7, 2, 1)
  brtsS   <- c(6, 5, 3, 2.5)
  tsplit  <- c(8)
  td      <- tsplit - 0.5
  pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
  pars2   <- c(100, 1, 1, tsplit, 0, 2)

  lik_P <- sls::loglik_slsP(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_P
  lik_Q <- sls::loglik_slsQ(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_Q

  testthat::expect_equal(
    lik_P, lik_Q
  )

  #test6
  lambdas <- c(0.2, 0.6)
  mus     <- c(0.1, 0.03)
  brtsM   <- c(10, 9.5, 8, 7, 2, 1)
  brtsS   <- c(6, 5, 3, 2.5)
  tsplit  <- c(8)
  td      <- tsplit - 0.5
  pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
  pars2   <- c(100, 1, 1, tsplit, 0, 2)

  lik_P <- sls::loglik_slsP(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_P
  lik_Q <- sls::loglik_slsQ(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_Q

  testthat::expect_true(
    abs(lik_P - lik_Q) < 1e-3
  )

  #test7
  lambdas <- c(0.2, 0.8)
  mus     <- c(0.1, 0.01)
  brtsM   <- c(10, 9.5, 8, 7, 2, 1)
  brtsS   <- c(6, 5, 3, 2.5)
  tsplit  <- c(8)
  td      <- tsplit - 0.5
  pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
  pars2   <- c(100, 1, 1, tsplit, 0, 2)

  lik_P <- sls::loglik_slsP(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_P
  lik_Q <- sls::loglik_slsQ(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS); lik_Q

  testthat::expect_true(
    abs(lik_P - lik_Q) < 1e-3
  )

})
