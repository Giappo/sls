context("conditional probabilities")

test_that("Analytical equivalence for conditional probability 1", {

  lambdas    <- c(0.5, 0.4)
  mus        <- c(0.4, 0.2)
  pars       <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  tp <- 0; tc <- -10; ts <- -6

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
  brtsM      <- c(10, 9, 1)
  brtsS      <- c(5, 3, 2)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  parsM <- c(lambdas[1], mus[1])
  parsS <- c(lambdas[2], mus[2])

  Pc1 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 1)
  Pc3 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 3)
  testthat::expect_true(Pc1 >= Pc3)

  #test2
  brtsM      <- c(13, 8, 3)
  brtsS      <- c(5)
  lambdas    <- c(0.5, 0.7)
  mus        <- c(0.4, 0.2)
  parsM <- c(lambdas[1], mus[1])
  parsS <- c(lambdas[2], mus[2])

  Pc1 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 1)
  Pc3 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 3)
  testthat::expect_true(Pc1 >= Pc3)

  #test3
  brtsM      <- c(10, 9, 6, 4, 1)
  brtsS      <- c(5, 3, 2.8)
  lambdas    <- c(0.4, 0.6)
  mus        <- c(0.3, 0.1)
  parsM <- c(lambdas[1], mus[1])
  parsS <- c(lambdas[2], mus[2])

  Pc1 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 1)
  Pc3 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 3)
  testthat::expect_true(Pc1 >= Pc3)
})

test_that("If lambda2=mu2=0 (inert subclade), Pc1 equal to Pc3 (PS = 1)", {

  #test1
  brtsM      <- c(10, 9, 1)
  brtsS      <- c(5, 3, 2)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  parsM <- c(lambdas[1], mus[1])
  parsS <- c(lambdas[2], mus[2])

  Pc1 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 1)
  Pc3 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 3)
  testthat::expect_true(Pc1 == Pc3)

  #test2
  brtsM      <- c(13, 8, 3)
  brtsS      <- c(5)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  parsM <- c(lambdas[1], mus[1])
  parsS <- c(lambdas[2], mus[2])

  Pc1 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 1)
  Pc3 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 3)
  testthat::expect_true(Pc1 == Pc3)

  #test3
  brtsM      <- c(10, 9, 6, 4, 1)
  brtsS      <- c(5, 3, 2.8)
  lambdas    <- c(0.4, 0)
  mus        <- c(0.3, 0)
  parsM <- c(lambdas[1], mus[1])
  parsS <- c(lambdas[2], mus[2])

  Pc1 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 1)
  Pc3 <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, parsM = parsM, parsS = parsS, cond = 3)
  testthat::expect_true(Pc1 == Pc3)
})

test_that("sls algorithm yields the same Pc1 provided by DDD", {

  testthat::skip('This test includes sls_sim, which I want to replace.')

  set.seed(s <- 2)
  parsM <- c(0.3, 0.1)
  parsS <- c(0.6, 0.08) #c(0.4, 0.2, 0.6, 0.1)

  cond <- 1
  age <- 10
  N0 <- 2
  t_d <- 4.8

  pars1 <- c(parsM[1], parsM[2], Inf, parsS[1], parsS[2], Inf, t_d)
  sim   <- sls::sls_sim(pars1 = pars1, age = age, N0 = N0, cond = cond)
  brts  <- sim$brts; brtsM <- brts[[1]]; brtsS <- brts[[2]]; brts

  tsplit <- min(abs(brtsM[abs(brtsM) > t_d]))
  missnumspec <- c(0, 0)
  res   <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))

  #pars2[3] is cond
  pars2_0 <- c(res, 1,    0, tsplit, 0, N0)
  pars2_1 <- c(res, 1, cond, tsplit, 0, N0)

  DDDloglik0 <- DDD::dd_KI_loglik(pars1 = pars1, pars2 = pars2_0,
                                  brtsM = abs(brtsM), brtsS = abs(brtsS), missnumspec = missnumspec)

  DDDloglik1 <- DDD::dd_KI_loglik(pars1 = pars1, pars2 = pars2_1,
                                  brtsM = abs(brtsM), brtsS = abs(brtsS), missnumspec = missnumspec)

  Pc <- sls::Pc_1shift(parsM = parsM, parsS = parsS, brtsM = brtsM, brtsS = c(t_d, brtsS),
                       cond = cond)

  test <- abs(DDDloglik1 - (DDDloglik0 - log(Pc)))

  testthat::expect_true(
    test < 1e-4 || (DDDloglik1 == -Inf && DDDloglik0 == -Inf)
  )

  if (packageVersion(pkg = 'DDD') >= 3.8)
  {
    cond <- 4
    pars2_0 <- c(res, 1,    0, tsplit, 0, N0)
    pars2_1 <- c(res, 1, cond, tsplit, 0, N0)

    DDDloglik0 <- DDD::dd_KI_loglik(pars1 = pars1, pars2 = pars2_0,
                                    brtsM = abs(brtsM), brtsS = abs(brtsS), missnumspec = missnumspec)

    DDDloglik1 <- DDD::dd_KI_loglik(pars1 = pars1, pars2 = pars2_1,
                                    brtsM = abs(brtsM), brtsS = abs(brtsS), missnumspec = missnumspec)

    Pc <- sls::Pc_1shift(parsM = parsM, parsS = parsS, brtsM = brtsM, brtsS = c(t_d, brtsS),
                         cond = cond)

    test <- abs(DDDloglik1 - (DDDloglik0 - log(Pc)))

    testthat::expect_true(
      test < 1e-4 || (DDDloglik1 == -Inf && DDDloglik0 == -Inf)
    )

  }

})
