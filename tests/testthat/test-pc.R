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
  brtsM      <- c(10, 9, 1)
  brtsS      <- c(3, 2)
  lambdas    <- c(0.5, 0)
  mus        <- c(0.4, 0)
  t_d <- 5
  pars1 <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, t_d)

  tsplit <- abs(max(brtsM[abs(brtsM) > t_d]))
  cond <- 1
  age <- 10
  soc <- 2
  missnumspec <- c(0,0)
  res   <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
  pars2 <- c(res, 1, cond, tsplit, 0 , soc)

  conditioning <- sls::Pc_1shift2(brtsM = brtsM, brtsS = brtsS, pars1 = pars1, pars2 = pars2)
  testthat::expect_true(conditioning[[1]] >= conditioning[[3]])

  #test2
  brtsM      <- c(13, 8, 3)
  brtsS      <- NULL
  lambdas    <- c(0.5, 0.7)
  mus        <- c(0.4, 0.2)
  t_d <- 5
  pars1 <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, t_d)

  tsplit <- abs(max(brtsM[abs(brtsM) > t_d]))
  cond <- 1
  age <- 10
  soc <- 2
  missnumspec <- c(0,0)
  res   <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
  pars2 <- c(res, 1, cond, tsplit, 0 , soc)

  conditioning <- sls::Pc_1shift2(brtsM = brtsM, brtsS = brtsS, pars1 = pars1, pars2 = pars2)
  testthat::expect_true(conditioning[[1]] >= conditioning[[3]])

  #test3
  brtsM      <- c(12, 9, 4, 2)
  brtsS      <- NULL
  lambdas    <- c(0.5, 0.4)
  mus        <- c(0.4, 0.3)
  t_d <- 5
  pars1 <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, t_d)

  tsplit <- abs(max(brtsM[abs(brtsM) > t_d]))
  cond <- 1
  age <- 10
  soc <- 2
  missnumspec <- c(0,0)
  res   <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
  pars2 <- c(res, 1, cond, tsplit, 0 , soc)

  conditioning <- sls::Pc_1shift2(brtsM = brtsM, brtsS = brtsS, pars1 = pars1, pars2 = pars2)
  testthat::expect_true(conditioning[[1]] >= conditioning[[3]])

  #test4
  brtsM      <- c(10, 9, 6, 4, 1)
  brtsS      <- c(3, 2.8)
  lambdas    <- c(0.4, 0.6)
  mus        <- c(0.3, 0.1)
  t_d <- 5
  pars1 <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, t_d)

  tsplit <- abs(max(brtsM[abs(brtsM) > t_d]))
  cond <- 1
  age <- 10
  soc <- 2
  missnumspec <- c(0,0)
  res   <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
  pars2 <- c(res, 1, cond, tsplit, 0 , soc)

  conditioning <- sls::Pc_1shift2(brtsM = brtsM, brtsS = brtsS, pars1 = pars1, pars2 = pars2)

  testthat::expect_true(conditioning[[1]] >= conditioning[[3]])
})

test_that("If lambda2=mu2=0 (inert subclade), Pc1 equal to Pc3 (PS = 1)", {

  brtsM      <- c(10, 9, 6, 4, 1)
  brtsS      <- c(3, 2.8)
  lambdas    <- c(0.4, 0)
  mus        <- c(0.3, 0)
  t_d <- 5
  pars1 <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, t_d)

  tsplit <- abs(max(brtsM[abs(brtsM) > t_d]))
  cond <- 1
  age <- 10
  soc <- 2
  missnumspec <- c(0,0)
  res   <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
  pars2 <- c(res, 1, cond, tsplit, 0 , soc)

  conditioning <- sls::Pc_1shift2(brtsM = brtsM, brtsS = brtsS, pars1 = pars1, pars2 = pars2)

  testthat::expect_true(conditioning[[1]] == conditioning[[3]])
})

test_that("sls algorithm yields the same Pc1 provided by DDD", {

  set.seed(s <- 2)
  simpars <- c(0.3, 0.1, 0.6, 0.08) #c(0.4, 0.2, 0.6, 0.1)

  cond <- 1
  age <- 10
  soc <- 2
  t_d <- 4.8

  pars1 <- c(simpars[1], simpars[2], Inf, simpars[3], simpars[4], Inf, t_d)
  sim   <- sls::sls_sim2(pars1 = pars1, age = age, soc = soc, cond = cond)
  brts  <- sim$brts; brtsM <- brts[[1]]; brtsS <- brts[[2]]; brts

  tsplit <- min(abs(brtsM[abs(brtsM) > t_d]))
  missnumspec <- c(0,0)
  res   <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))

  #pars2[3] is cond
  pars2_0 <- c(res, 1, 0, tsplit, 0 , soc)
  pars2_1 <- c(res, 1, cond, tsplit, 0 , soc)

  DDDloglik0 <- DDD::dd_KI_loglik(pars1 = pars1, pars2 = pars2_0,
                                  brtsM = abs(brtsM), brtsS = abs(brtsS), missnumspec = missnumspec)

  DDDloglik1 <- DDD::dd_KI_loglik(pars1 = pars1, pars2 = pars2_1,
                                  brtsM = abs(brtsM), brtsS = abs(brtsS), missnumspec = missnumspec)

  conditioning <- sls::Pc_1shift2(pars1 = pars1, pars2 = pars2_1, brtsM = brtsM, brtsS = brtsS)
  Pc <- conditioning[[cond]]; Pc

  test <- abs(DDDloglik1 - (DDDloglik0 - log(Pc)))

  testthat::expect_true(
    test < 1e-4 || (DDDloglik1 == -Inf && DDDloglik0 == -Inf)
  )

})
