context("likelihoods")

test_that( "test P and Q approach equivalence", {

  #test1
  brtsM      <- c(-12, -9)
  brtsS      <- c(-5)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P1 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q1 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_equal(
    lik_P1, lik_Q1
  )


  #test2
  brtsM      <- c(-12, -9, -4)
  brtsS      <- c(-5)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P2 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q2 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_equal(
    lik_P2, lik_Q2
  )

  #test3
  brtsM      <- c(-12, -9, -4)
  brtsS      <- c(-5)
  lambdas    <- c(0.5, 0.4)
  mus        <- c(0.4, 0.3)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P3 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q3 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_equal(
    lik_P3, lik_Q3
  )

  #test4
  brtsM      <- c(-12, -9, -4, -2)
  brtsS      <- c(-5)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P4 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q4 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_equal(
    lik_P4, lik_Q4
  )

  #test5
  brtsM      <- c(-12, -9, -4, -2)
  brtsS      <- c(-5)
  lambdas    <- c(0.5, 0.4)
  mus        <- c(0.4, 0.3)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P5 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q5 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_equal(
    lik_P5, lik_Q5
  )

  #test6
  brtsM      <- c(-12, -9, -4)
  brtsS      <- c(-5, -2)
  lambdas    <- c(0.5, 0.2)
  mus        <- c(0.4, 0.1)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P6 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q6 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_equal(
    lik_P6, lik_Q6
  )

  #test7
  brtsM      <- c(-12, -9, -4)
  brtsS      <- c(-5, -2)
  lambdas    <- c(0.5, 0.4)
  mus        <- c(0.4, 0.2)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P7 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q7 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_equal(
    lik_P7, lik_Q7
  )

  #test8
  brtsM      <- c(-12, -9, -4, -3)
  brtsS      <- c(-8, -2, -1)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P8 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q8 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_equal(
    lik_P8, lik_Q8
  )

  #test9
  brtsM      <- c(-12, -9, -4, -3)
  brtsS      <- c(-8, -2, -1)
  lambdas    <- c(0.5, 0.4)
  mus        <- c(0.4, 0.2)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P9 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q9 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)
  testthat::expect_equal(
    lik_P9, lik_Q9
  )

  #test10
  brtsM      <- c(-10, -10)
  brtsS      <- c(-4)
  lambdas    <- c(0.5, 0.4)
  mus        <- c(0.4, 0.2)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P10 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q10 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)

  testthat::expect_equal(
    lik_P10, lik_Q10
  )

  #test11
  brtsM      <- c(-9, -9)
  brtsS      <- c(-3)
  lambdas    <- c(0.4, 0.6)
  mus        <- c(0.3, 0.1)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P11 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q11 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)

  testthat::expect_equal(
    lik_P11, lik_Q11
  )

  #test12
  brtsM      <- c(-10, -10, -9, -6, -4, -1)
  brtsS      <- c(-5, -3)
  lambdas    <- c(0.4, 0.6)
  mus        <- c(0.3, 0.1)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P12 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q12 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)

  testthat::expect_true(
    abs(lik_P12 - lik_Q12) < 10^-3
  )

  #test13
  brtsM      <- c(-10, -10, -4)
  brtsS      <- c(-5, -2)
  lambdas    <- c(0.4, 0.1)
  mus        <- c(0.3, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P13 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q13 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars)

  testthat::expect_true(
    abs(lik_P13 - lik_Q13) < 10^-3
  )

  #test14
  brtsM      <- c(-10, -10, -4)
  brtsS      <- c(-5)
  lambdas    <- c(0.4, 0)
  mus        <- c(0.4, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  lik_P14 <- sls:::lik_shift_P(brtsM = brtsM, brtsS = brtsS, pars = pars)
  lik_Q14 <- loglik_Q <- sls:::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, pars = pars, lx = 80)

  testthat::expect_true(
    abs(lik_P14 - lik_Q14) < 10^-3
  )

})

test_that("test DDD_KI_loglik version", {

  #test1
  brtsM      <- c(-10, -10, -4)
  brtsS      <- c(-5)
  lambdas    <- c(0.4, 0.2)
  mus        <- c(0.3, 0.1)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  L0 <- sls:::lik_shift_DDD(pars = pars, brtsM = brtsM, brtsS = brtsS, cond = 0); L0
  L1 <- sls:::lik_shift_DDD(pars = pars, brtsM = brtsM, brtsS = brtsS, cond = 1); L1
  L2 <- sls:::lik_shift_DDD(pars = pars, brtsM = brtsM, brtsS = brtsS, cond = 2); L2
  L3 <- sls:::lik_shift_DDD(pars = pars, brtsM = brtsM, brtsS = brtsS, cond = 3); L3

  testthat::expect_true(L0 <= 0 & is.numeric(L0) & !is.infinite(L0))
  testthat::expect_true(L1 <= 0 & is.numeric(L1) & !is.infinite(L1))
  testthat::expect_true(L2 <= 0 & is.numeric(L2) & !is.infinite(L2))
  testthat::expect_true(L3 <= 0 & is.numeric(L3) & !is.infinite(L3))

  #test1
  brtsM      <- c(-10, -10, -2)
  brtsS      <- c(-4, -1)
  lambdas    <- c(0.4, 1e-5)
  mus        <- c(0.3, 0)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  L0 <- sls:::lik_shift_DDD(pars = pars, brtsM = brtsM, brtsS = brtsS, cond = 0); L0
  L1 <- sls:::lik_shift_DDD(pars = pars, brtsM = brtsM, brtsS = brtsS, cond = 1); L1
  L2 <- sls:::lik_shift_DDD(pars = pars, brtsM = brtsM, brtsS = brtsS, cond = 2); L2
  L3 <- sls:::lik_shift_DDD(pars = pars, brtsM = brtsM, brtsS = brtsS, cond = 3); L3

  testthat::expect_true(L0 <= 0 & is.numeric(L0) & !is.infinite(L0))
  testthat::expect_true(L1 <= 0 & is.numeric(L1) & !is.infinite(L1))
  testthat::expect_true(L2 <= 0 & is.numeric(L2) & !is.infinite(L2))
  testthat::expect_true(L3 <= 0 & is.numeric(L3) & !is.infinite(L3))
})

test_that("test sls_loglik_choosepar", {

  brtsM      <- c(-10, -10, -4)
  brtsS      <- c(-5, -3)
  lambdas    <- c(0.4, 0.2)
  mus        <- c(0.3, 0.1)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])

  #test1
  idparsfix  <- 1:2
  idparsopt  <- 3:4
  cond <- 2
  loglik_function <- sls::lik_shift_P
  test1a <- sls::sls_loglik_choosepar(loglik_function = loglik_function,
                                     brtsM = brtsM, brtsS = brtsS,
                                     idparsopt = idparsopt, idparsfix = idparsfix,
                                     trparsopt = pars[idparsopt],trparsfix = pars[idparsfix],
                                     cond = cond, pars.transform = 0)
  test1b <- loglik_function(brtsM = brtsM, brtsS = brtsS, pars = pars, cond = cond)
  testthat::expect_true(abs(test1a - test1b) < 1e-3)

  loglik_function <- sls:::lik_shift_DDD
  test2a <- sls::sls_loglik_choosepar(loglik_function = loglik_function,
                                      brtsM = brtsM, brtsS = brtsS,
                                      idparsopt = idparsopt, idparsfix = idparsfix,
                                      trparsopt = pars[idparsopt],trparsfix = pars[idparsfix],
                                      cond = cond, pars.transform = 0)
  test2b <- loglik_function(brtsM = brtsM, brtsS = brtsS, pars = pars, cond = cond)
  testthat::expect_true(abs(test2a - test2b) < 1e-3)

  #test2
  idparsfix  <- NULL
  idparsopt  <- 1:4
  cond <- 1
  loglik_function <- sls::lik_shift_P
  test1a <- sls::sls_loglik_choosepar(loglik_function = loglik_function,
                                      brtsM = brtsM, brtsS = brtsS,
                                      idparsopt = idparsopt, idparsfix = idparsfix,
                                      trparsopt = pars[idparsopt],trparsfix = pars[idparsfix],
                                      cond = cond, pars.transform = 0)
  test1b <- loglik_function(brtsM = brtsM, brtsS = brtsS, pars = pars, cond = cond)
  testthat::expect_true(abs(test1a - test1b) < 1e-3)

  loglik_function <- sls:::lik_shift_DDD
  test2a <- sls::sls_loglik_choosepar(loglik_function = loglik_function,
                                      brtsM = brtsM, brtsS = brtsS,
                                      idparsopt = idparsopt, idparsfix = idparsfix,
                                      trparsopt = pars[idparsopt],trparsfix = pars[idparsfix],
                                      cond = cond, pars.transform = 0)
  test2b <- loglik_function(brtsM = brtsM, brtsS = brtsS, pars = pars, cond = cond)
  testthat::expect_true(abs(test2a - test2b) < 1e-3)

})
