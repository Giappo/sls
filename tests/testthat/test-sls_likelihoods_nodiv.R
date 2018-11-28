context("likelihoods - no division")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("all the likelihoods with no division yield the same result", {

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

    res_1_1 <- fun1(
      pars_m = pars_m1,
      pars_s = pars_s1,
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_max = precision
    ); res_1_1
    res_2_1 <- fun2(
      pars_m = pars_m1,
      pars_s = pars_s1,
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_max = precision
    ); res_2_1

    res_1_2 <- res_2_2 <- 0
    if (ratios == TRUE) {
      pars_m2 <- pars_m / 2; pars_s2 <- pars_s * 3 / 4;

      res_1_2 <- fun1(
        pars_m = pars_m2,
        pars_s = pars_s2,
        brts_m = brts_m,
        brts_s = brts_s,
        cond = cond,
        n_max = precision
      ); res_1_2
      res_2_2 <- fun2(
        pars_m = pars_m2,
        pars_s = pars_s2,
        brts_m = brts_m,
        brts_s = brts_s,
        cond = cond,
        n_max = precision
      ); res_2_2
    }

    delta_1 <- res_1_1 - res_1_2; delta_1
    delta_2 <- res_2_1 - res_2_2; delta_2

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
    sls::loglik_bisse_shift,
    sls::loglik_ddd,
    sls::loglik_sls_p_nodiv,
    sls::loglik_sls_q_nodiv
  )
  threshold <- (!is_on_ci()) * 1e-2 +
               (is_on_ci())  * 1e-3

  cond <- 0; s <- 1
  for (s in 1:(2 + 4 * is_on_ci())) {
    set.seed(s)
    t_0s    <- c(4, 1.5)
    brts_m  <- c(
      t_0s[1],
      sort(runif(n = 20, min = 0.01, max = t_0s[1] - 0.01), decreasing = TRUE)
    )
    pars_m  <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    )
    brts_s  <- c(
      t_0s[2],
      sort(runif(n = 10, min = 0.01, max = t_0s[2] - 0.01), decreasing = TRUE)
    )
    pars_s  <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    ) * c(2, 0.5)
    cond   <- (cond == 0) * 1 + (cond == 1) * 0

    tests <- 0
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
        tests <- tests + 1
      }
    }
  }

  testthat::expect_equal(
    tests, length(models) * (length(models) - 1) / 2
  )
})
