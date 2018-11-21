context("likelihoods - bisse")

test_that( "test bisse and DDD logliks equivalence", {

  diff <- function(pars, brts) {

    bisse1 <- sls::loglik_bisse(pars, brts)
    bisse2 <- sls::loglik_bisse(pars / 2, brts)
    ddd1   <- DDD::bd_loglik(
      pars1 = c(pars, 0, 0),
      pars2 = c(0, 0, 1, 0, 2),
      brts = brts,
      missnumspec = 0
    )
    ddd2   <- DDD::bd_loglik(
      pars1 = c(pars / 2, 0, 0),
      pars2 = c(0, 0, 1, 0, 2),
      brts = brts,
      missnumspec = 0
    )

    Delta_ddd   <- ddd1 - ddd2
    Deltabisse <- bisse1 - bisse2

    diff <- abs(Deltabisse - Delta_ddd)

    return(diff)
  }

  #test1
  brts  <- c(10, 4, 2)
  pars  <- c(0.3, 0.1)

  testthat::expect_true(
    diff(pars = pars, brts = brts) <= 1e-5
  )

  #test2
  for (s in 1:20) {
    set.seed(s)
    brts  <- c(
      10,
      sort(runif(n = 30, min = 0.01, max = 10 - 0.01), decreasing = TRUE)
    )
    pars  <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    )

    testthat::expect_true(
      diff(pars = pars, brts = brts) <= 1e-5
    )
  }
})

test_that( "test bisse alternative functions", {

  diff4 <- function(pars, brts) {

    bisse_a1 <- sls::loglik_bisse(pars, brts)
    bisse_a2 <- sls::loglik_bisse(pars / 2, brts)
    bisse_b1 <- sls::loglik_bisse2(pars, brts)
    bisse_b2 <- sls::loglik_bisse2(pars / 2, brts)

    diff1 <- abs(bisse_a1 - bisse_b1)
    diff2 <- abs(bisse_a2 - bisse_b2)

    return(c(diff1, diff2))
  }

  #test1
  brts  <- c(10, 4, 2)
  pars  <- c(0.3, 0.1)

  testthat::expect_true(
    all(diff4(pars = pars, brts = brts) <= 1e-5)
  )

  #test2
  for (s in 1:20) {
    set.seed(s)
    brts  <- c(
      10,
      sort(runif(n = 30, min = 0.01, max = 10 - 0.01), decreasing = TRUE)
    )
    pars  <- c(
      x <- runif(n = 1, min = 0.1, max = 1),
      runif(n = 1, min = 0.05, max = x * 3 / 4)
    )

    testthat::expect_true(
      all(diff4(pars = pars, brts = brts) <= 1e-5)
    )
  }
})

test_that("test bisse alternative functions for the version with shift", {
  #  this is a test for the old wrong method

  diff5 <- function(
    pars_m,
    pars_s,
    brts_m,
    brts_s,
    cond,
    precision = 1e2
  ) {

    pars_m1 <- pars_m  ; pars_s1 <- pars_s;
    pars_m2 <- pars_m / 2; pars_s2 <- pars_s * 3 / 4;

    bisse_a1 <- sls::loglik_bisse_shift(
      pars_m = pars_m1,
      pars_s = pars_s1,
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_max = precision
    ); bisse_a1
    bisse_a2 <- sls::loglik_bisse_shift(
      pars_m = pars_m2,
      pars_s = pars_s2,
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_max = precision
    ); bisse_a2
    bisse_b1 <- sls::loglik_bisse_shift2(
      pars = pars_m1,
      brts = brts_m,
      t_d = brts_s[1]
    ) +
      sls::loglik_bisse2(
        pars = pars_s1,
        brts = brts_s,
        n_0 = 1
      ); bisse_b1
    bisse_b2 <- sls::loglik_bisse_shift2(
      pars = pars_m2,
      brts = brts_m,
      t_d = brts_s[1]
    ) +
      sls::loglik_bisse2(
        pars = pars_s2,
        brts = brts_s,
        n_0 = 1
      ); bisse_b2

    diff1 <- abs(bisse_a1 - bisse_b1)
    diff2 <- abs(bisse_a2 - bisse_b2)

    return(c(diff1, diff2))
  }

#test1
brts_m <- c(10, 4, 2)
pars_m <- c(0.3, 0.1)
brts_s <- c(3, 1, 0.5)
pars_s <- c(0.5, 0.05)
cond   <- 0

testthat::expect_true(
  all(diff5(
    pars_m = pars_m,
    brts_m = brts_m,
    pars_s = pars_s,
    brts_s = brts_s,
    cond = cond) <= 1e-5)
)

#test2
l_m <- 13 + (2 * (ribir::is_on_travis()))
age <- 10;
maxs <- 10 + (90 * (ribir::is_on_travis()))
res <- rep(NA, maxs)
test_threshold <- 1e-3
max_iterations <- 8 + (ribir::is_on_travis())
for (s in 1:maxs) {
  set.seed(s)

  l1  <- runif(n = 1, min = 0.1, max = 1)
  m1  <- runif(n = 1, min = 0.02, max = l1 * (3 / 4))
  l2  <- l1 * 2
  m2  <- m1 / 2
  pars_m   <- c(l1, m1)
  pars_s   <- c(l2, m2)
  brts_m   <- c(
    age,
    sort(runif(n = (l_m - 1), min = 0, max = age), decreasing = TRUE)
  )
  tsplit  <- sample(
    x = brts_m[-c(1:floor(l_m / 6), (l_m - floor(l_m / 6)):l_m)],
    size = 1
  )
  t_d      <- tsplit - 0.1
  brts_s  <- sort(
    runif(n = floor(l_m / 2), min = 0, max = t_d - 0.1),
    decreasing = TRUE
  )
  cond <- 0

  testthat::expect_true(
    all(diff5(
      pars_m = pars_m,
      brts_m = brts_m,
      pars_s = pars_s,
      brts_s = brts_s,
      cond = cond,
      precision = 1e2
    ) <= 1e-5)
  )
}

})
