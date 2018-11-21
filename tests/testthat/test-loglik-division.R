context("likelihoods - division")

test_that( "all the likelihoods with division yield the same result", {

  # testthat::skip('I skip it because it is slow. It works, though.') # nolint

  while (!require("ribir")) {
    devtools::install_github("richelbilderbeek/ribir")
  }

  diff <- function(
    pars_m,
    pars_s,
    brts_m,
    brts_s,
    cond,
    fun1,
    fun2,
    precision = 1e2,
    ratios = FALSE
  ) {
    pars_m1 <- pars_m  ; pars_s1 <- pars_s

    res1.1 <- fun1(
      pars_m = pars_m1,
      pars_s = pars_s1,
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_max = precision
    ); res1.1

    res2.1 <- fun2(
      pars_m = pars_m1,
      pars_s = pars_s1,
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_max = precision
    ); res2.1

    res1.2 <- res2.2 <- 0
    if (ratios == TRUE) {
      pars_m2 <- pars_m / 2; pars_s2 <- pars_s * 3 / 4;

      res1.2 <- fun1(
        pars_m = pars_m2,
        pars_s = pars_s2,
        brts_m = brts_m,
        brts_s = brts_s,
        cond = cond,
        n_max = precision
      ); res1.2

      res2.2 <- fun2(
        pars_m = pars_m2,
        pars_s = pars_s2,
        brts_m = brts_m,
        brts_s = brts_s,
        cond = cond,
        n_max = precision
      ); res2.2
    }

    delta_1 <- res1.1 - res1.2; delta_1
    delta_2 <- res2.1 - res2.2; delta_2

    diff <- abs(delta_1 - delta_2)

    return(diff)
  }
  test_diff <- function(
    pars_m,
    pars_s,
    brts_m,
    brts_s,
    cond,
    fun1,
    fun2,
    precision = 1e2,
    threshold = 1e-3
  ) {
    out <- 1; max_rep <- 10; rep <- 0
    while (out > threshold && rep <= max_rep) {
      precision <- precision * 2
      out <- diff(
        pars_m = pars_m,
        pars_s = pars_s,
        brts_m = brts_m,
        brts_s = brts_s,
        cond = cond,
        fun1 = fun1,
        fun2 = fun2,
        precision = precision
      )
      rep <- rep + 1
    }
    out
  }

  models <- c(
    sls::loglik_sls_p,
    sls::loglik_sls_q
  )
  threshold <- (!ribir::is_on_travis()) * 1e-2 +
               (ribir::is_on_travis())  * (1 / 2) * 1e-3

  cond <- sls_conds()[1]
  for (s in 1:(4 + 4 * ribir::is_on_travis())) {
    set.seed(s)
    t_0s    <- c(6, 2)
    brts_m  <- c(
      t_0s[1],
      sort(runif(n = 20, min = 0.01, max = t_0s[1] - 0.01), decreasing = TRUE)
    )
    pars_m  <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    )
    brts_s <- c(
      t_0s[2],
      sort(runif(n = 10, min = 0.01, max = t_0s[2] - 0.01), decreasing = TRUE)
    )
    pars_s <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    ) * c(2, 0.5)
    cond <- c(sls_conds(), sls_conds())[which(cond %in% sls_conds()) + 1]

    for (i in 1:(length(models) - 1)) {
      for (j in (i + 1):length(models)) {
        testthat::expect_true(
          test_diff(
            pars_m = pars_m,
            pars_s = pars_s,
            brts_m = brts_m,
            brts_s = brts_s,
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
