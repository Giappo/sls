context("likelihoods - bisse")

test_that( "test BISSE and DDD logliks equivalence", {

  diff <- function(pars, brts) {

    BISSE1 <- sls::loglik_bisse(pars, brts)
    BISSE2 <- sls::loglik_bisse(pars/2, brts)
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

test_that( "test BISSE alternative functions", {

  diff4 <- function(pars, brts) {

    BISSEA1 <- sls::loglik_bisse(pars, brts)
    BISSEA2 <- sls::loglik_bisse(pars/2, brts)
    BISSEB1 <- sls::loglik_bisse2(pars, brts)
    BISSEB2 <- sls::loglik_bisse2(pars/2, brts)

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

diff5 <- function(parsM, parsS, brtsM, brtsS, cond, precision = 1e2) {

  parsM1 <- parsM  ; parsS1 <- parsS;
  parsM2 <- parsM/2; parsS2 <- parsS * 3/4;

  BISSEA1 <- sls::loglik_bisse_shift(parsM = parsM1, parsS = parsS1,
                                    brtsM = brtsM, brtsS = brtsS, cond = cond,
                                    nmax = precision); BISSEA1
  BISSEA2 <- sls::loglik_bisse_shift(parsM = parsM2, parsS = parsS2,
                                    brtsM = brtsM, brtsS = brtsS, cond = cond,
                                    nmax = precision); BISSEA2
  BISSEB1 <- sls::loglik_bisse_shift2(pars = parsM1, brts = brtsM, td = brtsS[1]) +
             sls::loglik_bisse2(pars = parsS1, brts = brtsS, N0 = 1); BISSEB1
  BISSEB2 <- sls::loglik_bisse_shift2(pars = parsM2, brts = brtsM, td = brtsS[1]) +
             sls::loglik_bisse2(pars = parsS2, brts = brtsS, N0 = 1); BISSEB2


  diff1 <- abs(BISSEA1 - BISSEB1)
  diff2 <- abs(BISSEA2 - BISSEB2)

  return(c(diff1, diff2))
}

#test1
brtsM  <- c(10, 4, 2)
parsM  <- c(0.3, 0.1)
brtsS  <- c(3, 1, 0.5)
parsS  <- c(0.5, 0.05)
cond   <- 0

testthat::expect_true(
  all(diff5(parsM = parsM, brtsM = brtsM, parsS = parsS, brtsS = brtsS, cond = cond) <= 1e-5)
)

#test2
lM <- 13 + (2 * (ribir:::is_on_travis())); age <- 10;
maxs <- 10 + (90 * (ribir:::is_on_travis())); res <- rep(NA, maxs); test_threshold <- 1e-3; max_iterations <- 8 + (ribir:::is_on_travis())
for (s in 1:maxs)
{
  set.seed(s)
  # brts  <- c(10, sort(runif(n = 30, min = 0.01, max = 10 - 0.01), decreasing = TRUE))
  # pars  <- c(x <- runif(n = 1, min = 0.1, max = 1), runif(n = 1, min = 0.05, max = x*3/4))

  l1  <- runif(n = 1, min = 0.1 , max = 1)
  m1  <- runif(n = 1, min = 0.02, max = l1 * (3/4))
  l2  <- l1 * 2
  m2  <- m1 / 2
  parsM   <- c(l1, m1)
  parsS   <- c(l2, m2)
  brtsM   <- c(age, sort(runif(n = (lM - 1), min = 0, max = age), decreasing = TRUE))
  tsplit  <- sample(x = brtsM[-c(1:floor(lM/6), (lM - floor(lM/6)):lM)], size = 1)
  td      <- tsplit - 0.1
  brtsS   <- sort(runif(n = floor(lM/2), min = 0, max = td - 0.1), decreasing = TRUE)
  cond    <- 0

  testthat::expect_true(
    all(diff5(parsM = parsM, brtsM = brtsM, parsS = parsS, brtsS = brtsS, cond = cond, precision = 1e2) <= 1e-5)
  )
}

})
