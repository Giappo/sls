context("likelihoods - division")

test_that( "all the likelihoods with division yield the same result", {

  testthat::skip('I skip it because it is slow. It works, though.')

  while (!require("ribir")) {devtools::install_github("richelbilderbeek/ribir")}

  diff <- function(
    parsM, parsS, brtsM, brtsS, cond,
    fun1, fun2,
    precision = 1e2,
    ratios = FALSE
  )
  {

    parsM1 <- parsM  ; parsS1 <- parsS

    res1.1 <- fun1(parsM = parsM1,
                   parsS = parsS1,
                   brtsM = brtsM,
                   brtsS = brtsS,
                   cond = cond,
                   nmax = precision
    ); res1.1

    res2.1 <- fun2(
      parsM = parsM1,
      parsS = parsS1,
      brtsM = brtsM,
      brtsS = brtsS,
      cond = cond,
      nmax = precision
    ); res2.1

    res1.2 <- res2.2 <- 0
    if (ratios == TRUE)
    {
      parsM2 <- parsM/2; parsS2 <- parsS * 3/4;

      res1.2 <- fun1(parsM = parsM2,
                     parsS = parsS2,
                     brtsM = brtsM,
                     brtsS = brtsS,
                     cond = cond,
                     nmax = precision
      ); res1.2

      res2.2 <- fun2(
        parsM = parsM2,
        parsS = parsS2,
        brtsM = brtsM,
        brtsS = brtsS,
        cond = cond,
        nmax = precision
      ); res2.2
    }

    Delta1 <- res1.1 - res1.2; Delta1
    Delta2 <- res2.1 - res2.2; Delta2

    diff <- abs(Delta1 - Delta2)

    return(diff)
  }
  test.diff <- function(
    parsM, parsS, brtsM, brtsS, cond,
    fun1, fun2,
    precision = 1e2, threshold = 1e-3
  )
  {
    out <- 1; max_rep <- 10; rep <- 0
    while (out > threshold && rep <= max_rep)
    {
      precision <- precision * 2
      out <- diff(parsM = parsM, parsS = parsS, brtsM = brtsM, brtsS = brtsS, cond = cond,
                  fun1 = fun1, fun2 = fun2,
                  precision = precision)
      rep <- rep + 1
    }
    out
  }

  models <- c(sls::loglik_slsP,
              sls::loglik_slsQ
  )
  threshold <- (!ribir::is_on_travis()) * 10^-2 +
               (ribir::is_on_travis())  * (1/2) * 10^-3

  cond <- 0
  for (s in 1:(4 + 4 * ribir::is_on_travis()))
  {
    set.seed(s)
    t0s    <- c(6, 2)
    brtsM  <- c(t0s[1], sort(runif(n = 20, min = 0.01, max = t0s[1] - 0.01), decreasing = TRUE))
    parsM  <- c(x <- runif(n = 1, min = 0.1, max = 1), runif(n = 1, min = 0.05, max = x * 3/4))
    brtsS  <- c(t0s[2], sort(runif(n = 10, min = 0.01, max = t0s[2] - 0.01), decreasing = TRUE))
    parsS  <- c(x <- runif(n = 1, min = 0.1, max = 1), runif(n = 1, min = 0.05, max = x * 3/4)) * c(2, 0.5)
    cond   <- (cond == 0) * 1 + (cond == 1) * 0

    for (i in 1:(length(models) - 1))
    {
      for (j in (i + 1):length(models))
      {
        testthat::expect_true(
          test.diff(
            parsM = parsM,
            parsS = parsS,
            brtsM = brtsM,
            brtsS = brtsS,
            cond = cond,
            fun1 = models[[i]],
            fun2 = models[[j]],
            threshold = threshold
          ) <= threshold
        )
      }
    }
  }

})
