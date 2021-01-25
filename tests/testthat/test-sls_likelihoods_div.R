context("likelihoods - division")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}
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

test_that("all the likelihoods with division yield the same result", {

  models <- c(
    sls::loglik_sls_p,
    sls::loglik_sls_q
  )
  threshold <- (!is_on_ci()) * 1e-2 +
               (is_on_ci())  * (1 / 2) * 1e-3

  conds <- sls::sls_conds()
  cond <- conds[1]
  t_0s <- c(5, 2)
  n_m <- 15
  n_s <- 7
  for (seed in 1:(4 + 6 * is_on_ci())) {
    set.seed(seed)
    brts_m  <- c(
      t_0s[1],
      sort(runif(n = n_m, min = 0.01, max = t_0s[1] - 0.01), decreasing = TRUE)
    )
    pars_m  <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    )
    brts_s <- c(
      t_0s[2],
      sort(runif(n = n_s, min = 0.01, max = t_0s[2] - 0.01), decreasing = TRUE)
    )
    pars_s <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    ) * c(2, 0.5)
    pars <- c(pars_m, pars_s)
    brts <- list(brts_m, brts_s)
    cond <- c(conds, conds)[which(cond %in% conds) + 1]

    for (i in 1:(length(models) - 1)) {
      for (j in (i + 1):length(models)) {
        fun1 <- models[[i]]
        fun2 <- models[[j]]
        diff <- 1; max_rep <- 20; rep <- 0; precision <- 1e2
        while (diff > threshold && rep <= max_rep) {
          precision <- precision * 2
          res_1_1 <- fun1(
            pars = pars,
            brts = brts,
            cond = cond,
            n_max = precision
          ); res_1_1
          res_2_1 <- fun2(
            pars = pars,
            brts = brts,
            cond = cond,
            n_max = precision
          ); res_2_1

          diff <- abs(res_1_1 - res_2_1)
          rep <- rep + 1
        }
        testthat::expect_lt(diff, threshold)
      }
    }
  }
})

test_that("faster likelihood gives the same results", {

  models <- c(
    sls::loglik_sls_p,
    sls::loglik_sls_p2
  )
  threshold <- (!is_on_ci()) * 1e-2 +
    (is_on_ci())  * (1 / 2) * 1e-3

  conds <- sls::sls_conds()
  cond <- conds[1]
  t_0s <- c(5, 2)
  n_m <- 15
  n_s <- 7
  for (seed in 1:(4 + 6 * is_on_ci())) {
    set.seed(seed)
    brts_m  <- c(
      t_0s[1],
      sort(runif(n = n_m, min = 0.01, max = t_0s[1] - 0.01), decreasing = TRUE)
    )
    pars_m  <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    )
    brts_s <- c(
      t_0s[2],
      sort(runif(n = n_s, min = 0.01, max = t_0s[2] - 0.01), decreasing = TRUE)
    )
    pars_s <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    ) * c(2, 0.5)
    pars <- c(pars_m, pars_s)
    brts <- list(brts_m, brts_s)
    cond <- c(conds, conds)[which(cond %in% conds) + 1]

    for (i in 1:(length(models) - 1)) {
      for (j in (i + 1):length(models)) {
        fun1 <- models[[i]]
        fun2 <- models[[j]]
        diff <- 1; max_rep <- 20; rep <- 0; precision <- 1e2
        while (diff > threshold && rep <= max_rep) {
          precision <- precision * 2
          res_1_1 <- fun1(
            pars = pars,
            brts = brts,
            cond = cond,
            n_max = precision
          ); res_1_1
          res_2_1 <- fun2(
            pars = pars,
            brts = brts,
            cond = cond,
            n_max = precision
          ); res_2_1

          diff <- abs(res_1_1 - res_2_1)
          rep <- rep + 1
        }
        testthat::expect_lt(diff, threshold)
      }
    }
  }
})

test_that("div and nodiv yield the same values for mu = 0", {
  n_0 <- 2
  n_max <- 1e2
  cond <- 2
  starting_seed <- 1000
  t_0s <- c(5, 2)
  for (seed in starting_seed:(starting_seed + 100 + 400 * is_on_ci())) {
    set.seed(seed)
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
    pars_m[2] <- pars_s[2] <- 0
    pars <- c(pars_m, pars_s)
    brts <- list(brts_m, brts_s)
    cond <- 2 * (cond == 3) + 3 * (cond == 2)
    testthat::expect_equal(
      sls::loglik_sls_p(
        pars = pars,
        brts = brts,
        cond = cond,
        n_0 = n_0,
        n_max = n_max
      ),
      sls::loglik_sls_p_nodiv(
        pars = pars,
        brts = brts,
        cond = cond,
        n_0 = n_0,
        n_max = n_max
      )
    )
  }
})
