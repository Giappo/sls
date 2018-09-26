context("likelihoods2")

test_that( "test P and Q approach equivalence", {

  #define difference
  likelihood_diff <- function(fun1 = sls::lik_shift_P2_noratio,
                              fun2 = DDD::dd_KI_loglik,
                              pars1,
                              pars2,
                              brtsM,
                              brtsS,
                              missnumspec) {

    lik1a <- fun1(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS,
                                       missnumspec = missnumspec); lik1a
    lik2a <- fun2(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS,
                               missnumspec = missnumspec); lik2a

    pars1b <- c(pars1[1:6]/2, pars1[7])

    lik1b <- fun1(pars1 = pars1b, pars2 = pars2 , brtsM = brtsM, brtsS = brtsS,
                                       missnumspec = missnumspec); lik1b
    lik2b <- fun2(pars1 = pars1b, pars2 = pars2, brtsM = brtsM, brtsS = brtsS,
                               missnumspec = missnumspec); lik2b

    lik1 <- lik1a - lik1b
    lik2 <- lik2a - lik2b

    return(lik1 - lik2)
  }

  #test1
  if (1) {
    lambdas <- c(0.3, 0.5)
    mus     <- c(0.1, 0.1)
    brtsM   <- c(10, 8, 7, 2)
    brtsS   <- c(3)
    tsplit  <- c(7)
    td      <- c(6.5)
    cond    <- 0
    pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
    pars2   <- c(800, 1, cond, tsplit, 0, 2)
    missnumspec <- 0
  }

  diff <- likelihood_diff(fun1 = sls::lik_shift_P2_noratio,
                          fun2 = DDD::dd_KI_loglik,
                          pars1 = pars1,
                          pars2 = pars2,
                          brtsM = brtsM,
                          brtsS = brtsS,
                          missnumspec = missnumspec); diff

  testthat::expect_true(
    abs(diff) < 1e-3
  )

  #test2
  if (1) {
    lambdas <- c(0.4, 0.6)
    mus     <- c(0.2, 0.1)
    brtsM   <- c(10, 8, 7, 5, 2, 1)
    brtsS   <- c(3, 2.8, 1.5)
    tsplit  <- c(5)
    td      <- c(4.8)
    cond    <- 0
    pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
    pars2   <- c(800, 1, cond, tsplit, 0, 2)
    missnumspec <- 0
  }

  diff <- likelihood_diff(fun1 = sls::lik_shift_P2_noratio,
                          fun2 = DDD::dd_KI_loglik,
                          pars1 = pars1,
                          pars2 = pars2,
                          brtsM = brtsM,
                          brtsS = brtsS,
                          missnumspec = missnumspec); diff

  testthat::expect_true(
    abs(diff) < 1e-3
  )

  #test3
  maxs <- 20; res <- rep(NA, maxs)
  for (s in 1:maxs) {
    set.seed(s)
    if (1) {
      l1 <- runif(n = 1, min = 0.1 , max = 1)
      m1 <- runif(n = 1, min = 0.02, max = l1*(3/4))
      l2 <- l1 * 2
      m2 <- m1 / 2
      lambdas <- c(l1, l2)
      mus     <- c(m1, m2)
      lM      <- 20
      brtsM   <- c(10, sort(runif(n = (lM - 1), min = 0, max = 10), decreasing = TRUE))
      tsplit  <- sample(x = brtsM[-c(1:3, 22:25)], size = 1)
      td      <- tsplit - 0.1
      brtsS   <- sort(runif(n = floor(lM)/2, min = 0, max = td - 0.1), decreasing = TRUE)
      cond    <- 0
      pars1   <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, td)
      pars2   <- c(500, 1, cond, tsplit, 0, 2)
      missnumspec <- 0
    }

    diff <- likelihood_diff(fun1 = sls::lik_shift_P2_noratio,
                            fun2 = DDD::dd_KI_loglik,
                            pars1 = pars1,
                            pars2 = pars2,
                            brtsM = brtsM,
                            brtsS = brtsS,
                            missnumspec = missnumspec); diff

    res[s] <- diff
  }

  testthat::expect_true(
    abs(res) < 1e-3 | is.infinite(abs(res))
  )

})
