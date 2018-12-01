context("likelihoods - division")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("all the likelihoods with division yield the same result", {

  # testthat::skip('I skip it because it is slow. It works, though.') # nolint

  diff <- function(
    pars,
    brts,
    cond,
    fun1,
    fun2,
    precision = 1e2,
    ratios = FALSE
  ) {
    pars_m <- pars[1:2]
    pars_s <- pars[3:4]
    brts_m <- brts[[1]]
    brts_s <- brts[[2]]
    pars_m1 <- pars_m  ; pars_s1 <- pars_s
    pars_1 <- c(pars_m1, pars_s1)
    brts <- list(brts_m, brts_s)

    res_1_1 <- fun1(
      pars = pars_1,
      brts = brts,
      cond = cond,
      n_max = precision
    ); res_1_1

    res_2_1 <- fun2(
      pars = pars_1,
      brts = brts,
      cond = cond,
      n_max = precision
    ); res_2_1

    res_1_2 <- res_2_2 <- 0
    if (ratios == TRUE) {
      pars_m2 <- pars_m / 2; pars_s2 <- pars_s * 3 / 4;
      pars_2 <- c(pars_m2, pars_s2)

      res_1_2 <- fun1(
        pars = pars_2,
        brts = brts,
        cond = cond,
        n_max = precision
      ); res_1_2

      res_2_2 <- fun2(
        pars = pars_2,
        brts = brts,
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
    pars,
    brts,
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
        pars = pars,
        brts = brts,
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
  threshold <- (!is_on_ci()) * 1e-2 +
               (is_on_ci())  * (1 / 2) * 1e-3

  cond <- sls_conds()[1]
  for (s in 1:(4 + 4 * is_on_ci())) {
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
    pars <- c(pars_m, pars_s)
    brts <- list(brts_m, brts_s)
    cond <- c(sls_conds(), sls_conds())[which(cond %in% sls_conds()) + 1]

    for (i in 1:(length(models) - 1)) {
      for (j in (i + 1):length(models)) {
        testthat::expect_true(
          test_diff(
            pars = pars,
            brts = brts,
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
