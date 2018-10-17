context("BISSE likelihoods")

test_that( "test BISSE and DDD logliks equivalence", {

  diff <- function(pars, brts) {

    BISSE1 <- sls::BISSE_loglik(pars, brts)
    BISSE2 <- sls::BISSE_loglik(pars/2, brts)
    DDD1   <- DDD::bd_loglik(pars1 = c(pars, 0, 0)  , pars2 = c(0, 0, 1, 0, 2), brts = brts, missnumspec = 0)
    DDD2   <- DDD::bd_loglik(pars1 = c(pars/2, 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brts, missnumspec = 0)

    DeltaDDD   <- DDD1 - DDD2
    DeltaBISSE <- BISSE1 - BISSE2

    diff <- abs(DeltaBISSE - DeltaDDD)

    return(diff)
  }

  #test1
  brts  <- c(10, 4, 2)
  pars  <- c(0.3, 0.1)

  testthat::expect_true(
    diff(pars = pars, brts = brts) <= 1e-5
  )

  #test2
  for (s in 1:20)
  {
    set.seed(s)
    brts  <- c(10, sort(runif(n = 30, min = 0.01, max = 10 - 0.01), decreasing = TRUE))
    pars  <- c(x <- runif(n = 1, min = 0.1, max = 1), runif(n = 1, min = 0.05, max = x*3/4))

    testthat::expect_true(
      diff(pars = pars, brts = brts) <= 1e-5
    )
  }

})

test_that( "test BISSE and DDD logliks equivalence in case of one shift", {

  diff2 <- function(parsM, parsS, brtsM, brtsS, res = 1000) {

    parsM1 <- parsM  ; parsS1 <- parsS;
    parsM2 <- parsM/2; parsS2 <- parsS * 3/4;

    BISSE1 <- sls::BISSE_loglik_shift(parsM = parsM1, parsS = parsS1,
                                      brtsM = brtsM, brtsS = brtsS); BISSE1
    BISSE2 <- sls::BISSE_loglik_shift(parsM = parsM2, parsS = parsS2,
                                      brtsM = brtsM, brtsS = brtsS); BISSE2

    pars1.1 <- c(parsM1[1], parsM1[2], Inf, parsS1[1], parsS1[2], Inf, brtsS[1])
    pars2.1 <- c(res, 1, 0, min(abs(brtsM[abs(brtsM) > pars1.1[7]])), 0, 2)
    DDD1   <- DDD::dd_KI_loglik(pars1 = pars1.1,
                                pars2 = pars2.1,
                                brtsM = brtsM,
                                brtsS = brtsS[-1],
                                missnumspec = c(0,0)
    ); DDD1

    pars1.2 <- c(parsM2[1], parsM2[2], Inf, parsS2[1], parsS2[2], Inf, brtsS[1])
    pars2.2 <- c(res, 1, 0, min(abs(brtsM[abs(brtsM) > pars1.2[7]])), 0, 2)
    DDD2   <- DDD::dd_KI_loglik(pars1 = pars1.2,
                                pars2 = pars2.2,
                                brtsM = brtsM,
                                brtsS = brtsS[-1],
                                missnumspec = c(0,0)
    ); DDD2

    DeltaDDD   <- DDD1   - DDD2
    DeltaBISSE <- BISSE1 - BISSE2

    diff <- abs(DeltaBISSE - DeltaDDD)

    return(diff)
  }

  diff3 <- function(parsM, parsS, brtsM, brtsS) {

    parsM1 <- parsM  ; parsS1 <- parsS;
    parsM2 <- parsM/2; parsS2 <- parsS * 3/4;

    BISSE1 <- sls::BISSE_loglik_shift(parsM = parsM1, parsS = parsS1,
                                      brtsM = brtsM, brtsS = brtsS); BISSE1
    BISSE2 <- sls::BISSE_loglik_shift(parsM = parsM2, parsS = parsS2,
                                      brtsM = brtsM, brtsS = brtsS); BISSE2

    res <- 200

    pars1.1 <- c(parsM1[1], parsM1[2], Inf, parsS1[1], parsS1[2], Inf, brtsS[1])
    pars2.1 <- c(res, 1, 0, min(abs(brtsM[abs(brtsM) > pars1.1[7]])), 0, 2)
    nodiv1 <- sls::loglik_slsP_nodivision(pars1 = pars1.1,
                                pars2 = pars2.1,
                                brtsM = brtsM,
                                brtsS = brtsS[-1],
                                missnumspec = c(0,0)
    ); nodiv1

    pars1.2 <- c(parsM2[1], parsM2[2], Inf, parsS2[1], parsS2[2], Inf, brtsS[1])
    pars2.2 <- c(res, 1, 0, min(abs(brtsM[abs(brtsM) > pars1.2[7]])), 0, 2)
    nodiv2 <- sls::loglik_slsP_nodivision(pars1 = pars1.2,
                                pars2 = pars2.2,
                                brtsM = brtsM,
                                brtsS = brtsS[-1],
                                missnumspec = c(0,0)
    ); nodiv2

    Deltanodiv <- nodiv1 - nodiv2; Deltanodiv
    DeltaBISSE <- BISSE1 - BISSE2; DeltaBISSE

    diff <- abs(DeltaBISSE - Deltanodiv)

    return(diff)
  }

  #test1
  brtsM <- c(10, 8, 7, 4, 2)
  brtsS <- c(3, 1, 0.5)
  parsM <- c(0.3, 0.1)
  parsS <- c(0.6, 0.05)

  testthat::expect_true(
    diff2(parsM = parsM, parsS = parsS, brtsM = brtsM, brtsS = brtsS) <= 1e-3
  )
  testthat::expect_true(
    diff3(parsM = parsM, parsS = parsS, brtsM = brtsM, brtsS = brtsS) <= 1e-5
  )

  #test2
  for (s in 1:20)
  {
    set.seed(s)
    t0s    <- c(10, 4)
    brtsM  <- c(t0s[1], sort(runif(n = 20, min = 0.01, max = t0s[1] - 0.01), decreasing = TRUE))
    parsM  <- c(x <- runif(n = 1, min = 0.1, max = 1), runif(n = 1, min = 0.05, max = x * 3/4))
    brtsS  <- c(t0s[2], sort(runif(n = 10, min = 0.01, max = t0s[2] - 0.01), decreasing = TRUE))
    parsS  <- c(x <- runif(n = 1, min = 0.1, max = 1), runif(n = 1, min = 0.05, max = x * 3/4)) * c(2, 0.5)

    # if (s <= 5) #this is not precise enough!
    # {
    #   delta2 <- diff2(parsM = parsM, parsS = parsS, brtsM = brtsM, brtsS = brtsS)
    #   print(delta2)
    #   testthat::expect_true(
    #     delta2 <= 1e-3
    #   )
    # }

    delta3 <- diff3(parsM = parsM, parsS = parsS, brtsM = brtsM, brtsS = brtsS)
    # print(delta3)
    testthat::expect_true(
      delta3 <= 1e-6
    )

  }

})

test_that( "test BISSE alternative functions", {

  diff4 <- function(pars, brts) {

    BISSEA1 <- sls::BISSE_loglik(pars, brts)
    BISSEA2 <- sls::BISSE_loglik(pars/2, brts)
    BISSEB1 <- sls::BISSE_loglik2(pars, brts)
    BISSEB2 <- sls::BISSE_loglik2(pars/2, brts)

    diff1 <- abs(BISSEA1 - BISSEB1)
    diff2 <- abs(BISSEA2 - BISSEB2)

    return(c(diff1, diff2))
  }

  #test1
  brts  <- c(10, 4, 2)
  pars  <- c(0.3, 0.1)

  testthat::expect_true(
    all(diff4(pars = pars, brts = brts) <= 1e-5)
  )

  #test2
  for (s in 1:20)
  {
    set.seed(s)
    brts  <- c(10, sort(runif(n = 30, min = 0.01, max = 10 - 0.01), decreasing = TRUE))
    pars  <- c(x <- runif(n = 1, min = 0.1, max = 1), runif(n = 1, min = 0.05, max = x*3/4))

    testthat::expect_true(
      all(diff4(pars = pars, brts = brts) <= 1e-5)
    )
  }

})

test_that( "test BISSE alternative functions for the version with shift (the old wrong one)", {

diff5 <- function(parsM, parsS, brtsM, brtsS) {

  parsM1 <- parsM  ; parsS1 <- parsS;
  parsM2 <- parsM/2; parsS2 <- parsS * 3/4;

  BISSEA1 <- sls::BISSE_loglik_shift(parsM = parsM1, parsS = parsS1,
                                    brtsM = brtsM, brtsS = brtsS); BISSEA1
  BISSEA2 <- sls::BISSE_loglik_shift(parsM = parsM2, parsS = parsS2,
                                    brtsM = brtsM, brtsS = brtsS); BISSEA2
  BISSEB1 <- sls::BISSE_loglik_shift2(pars = parsM1, brts = brtsM, td = brtsS[1]) +
             sls::BISSE_loglik2(pars = parsS1, brts = brtsS, N0 = 1); BISSEB1
  BISSEB2 <- sls::BISSE_loglik_shift2(pars = parsM2, brts = brtsM, td = brtsS[1]) +
             sls::BISSE_loglik2(pars = parsS2, brts = brtsS, N0 = 1); BISSEB2


  diff1 <- abs(BISSEA1 - BISSEB1)
  diff2 <- abs(BISSEA2 - BISSEB2)

  return(c(diff1, diff2))
}

#test1
brts  <- c(10, 4, 2)
pars  <- c(0.3, 0.1)

testthat::expect_true(
  all(diff5(parsM = parsM, brtsM = brtsM, parsS = parsS, brtsS = brtsS) <= 1e-5)
)

#test2
for (s in 1:20)
{
  set.seed(s)
  brts  <- c(10, sort(runif(n = 30, min = 0.01, max = 10 - 0.01), decreasing = TRUE))
  pars  <- c(x <- runif(n = 1, min = 0.1, max = 1), runif(n = 1, min = 0.05, max = x*3/4))

  testthat::expect_true(
    all(diff5(parsM = parsM, brtsM = brtsM, parsS = parsS, brtsS = brtsS) <= 1e-5)
  )
}

})
