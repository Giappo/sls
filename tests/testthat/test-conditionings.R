context("conditional probabilities")

test_that("Analytical equivalence for conditional probability 1", {

  lambdas <- c(0.5, 0.4)
  mus     <- c(0.4, 0.2)
  pars    <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  t_p <- 0; t_c <- -10; t_d <- -6

  ns <- 2:1e5
  p_ns_cs <- sls::pn(n = ns, t = t_d - t_c, lambda = lambdas[1], mu = mus[1])
  p_0_sp <- sls::pn(n = 0, t = t_p - t_d, lambda = lambdas[1], mu = mus[1])
  test1 <- sum(
    p_ns_cs * (1 - (p_0_sp) ^ (ns - 1) )
  )

  p_1_cs <- sls::pn(n = 1, t = t_d - t_c, lambda = lambdas[1], mu = mus[1])
  p_0_sp <- sls::pn(n = 0, t = t_p - t_d, lambda = lambdas[1], mu = mus[1])
  u_cs <- sls::ut(t = t_d - t_c, lambda = lambdas[1], mu = mus[1])
  p_0_sp <- sls::pn(n = 0, t = t_p - t_d, lambda = lambdas[1], mu = mus[1])
  test2 <- p_1_cs *
    (1 - p_0_sp) *
    u_cs *
    (1 - u_cs) ^ -1 *
    (1 - u_cs * p_0_sp) ^ -1

  testthat::expect_equal(
    test1, test2
  )

})

test_that("pc1 never smaller than pc3 (pc3 is a stricter condition)", {

  #test1
  brts_m  <- c(10, 9, 1)
  brts_s  <- c(5, 3, 2)
  lambdas <- c(0.5, 0)
  mus     <- c(0.4, 0)
  pars_m <- c(lambdas[1], mus[1])
  pars_s <- c(lambdas[2], mus[2])

  pc1 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 1
  )
  pc3 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 3
  )
  testthat::expect_true(pc1 >= pc3)

  #test2
  brts_m  <- c(13, 8, 3)
  brts_s  <- c(5)
  lambdas <- c(0.5, 0.7)
  mus     <- c(0.4, 0.2)
  pars_m <- c(lambdas[1], mus[1])
  pars_s <- c(lambdas[2], mus[2])

  pc1 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 1
  )
  pc3 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 3
  )
  testthat::expect_true(pc1 >= pc3)

  #test3
  brts_m  <- c(10, 9, 6, 4, 1)
  brts_s  <- c(5, 3, 2.8)
  lambdas <- c(0.4, 0.6)
  mus     <- c(0.3, 0.1)
  pars_m <- c(lambdas[1], mus[1])
  pars_s <- c(lambdas[2], mus[2])

  pc1 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 1
  )
  pc3 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 3
  )
  testthat::expect_true(pc1 >= pc3)
})

test_that("If lambda2=mu2=0 (inert subclade), pc1 equal to pc3 (PS = 1)", {

  #test1
  brts_m  <- c(10, 9, 1)
  brts_s  <- c(5, 3, 2)
  lambdas <- c(0.5, 0)
  mus     <- c(0.4, 0)
  pars_m <- c(lambdas[1], mus[1])
  pars_s <- c(lambdas[2], mus[2])

  pc1 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 1
  )
  pc3 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 3
  )
  testthat::expect_true(pc1 == pc3)

  #test2
  brts_m  <- c(13, 8, 3)
  brts_s  <- c(5)
  lambdas <- c(0.5, 0)
  mus     <- c(0.4, 0)
  pars_m <- c(lambdas[1], mus[1])
  pars_s <- c(lambdas[2], mus[2])

  pc1 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 1
  )
  pc3 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 3
  )
  testthat::expect_true(pc1 == pc3)

  #test3
  brts_m  <- c(10, 9, 6, 4, 1)
  brts_s  <- c(5, 3, 2.8)
  lambdas <- c(0.4, 0)
  mus     <- c(0.3, 0)
  pars_m <- c(lambdas[1], mus[1])
  pars_s <- c(lambdas[2], mus[2])

  pc1 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 1
  )
  pc3 <- sls::pc_1shift(
    brts_m = brts_m,
    brts_s = brts_s,
    pars_m = pars_m,
    pars_s = pars_s,
    cond = 3
  )
  testthat::expect_true(pc1 == pc3)
})

test_that("sls algorithm yields the same pc1 provided by DDD", {

  diff_pc_vs_ddd <- function(
    pars_m,
    pars_s,
    cond = 1,
    age,
    n_0,
    t_d,
    seed = 2,
    precision_threshold = 1e-4
  ) {

    lambdas <- c(pars_m[1], pars_s[1])
    mus <- c(pars_m[2], pars_s[2])

    set.seed(seed)
    sim <- sls::sls_sim(
      lambdas = lambdas,
      mus = mus,
      cond = cond
    )
    brts <- sim$brts; brts_m <- brts[[1]]; brts_s <- brts[[2]]; brts
    tsplit <- min(abs(brts_m[abs(brts_m) > t_d]))

    missnumspec <- c(0, 0)
    max_res <- 1200
    res <- min(10 * (1 + length(c(brts_m, brts_s)) + sum(missnumspec)), max_res)
    out <- 1
    while (out > precision_threshold && res <= max_res) {

      # pars2[3] is cond
      pars1 <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, t_d)
      pars2_0 <- c(res, 1,    0, tsplit, 0, n_0)
      pars2_1 <- c(res, 1, cond, tsplit, 0, n_0)

      ddd_loglik0 <- DDD::dd_KI_loglik(
        pars1 = pars1,
        pars2 = pars2_0,
        brtsM = brts_m,
        brtsS = brts_s,
        missnumspec = missnumspec
      )

      ddd_loglik1 <- DDD::dd_KI_loglik(
        pars1 = pars1,
        pars2 = pars2_1,
        brtsM = brts_m,
        brtsS = brts_s,
        missnumspec = missnumspec
      )

      pc <- sls::pc_1shift(
        pars_m = pars_m,
        pars_s = pars_s,
        brts_m = brts_m,
        brts_s = c(t_d, brts_s),
        cond = cond
      )

      ddd_pc <- (ddd_loglik0 - ddd_loglik1)
      out <- abs(
        (exp(ddd_pc) - exp(log(pc))) / exp(ddd_pc)
      )
      res <- 2 * res
    }
    out
  }

  pars_m <- c(0.3, 0.1)
  pars_s <- c(0.6, 0.08)
  age <- 10
  n_0 <- 2
  t_d <- 4.8
  precision_threshold <- 1e-3

  test <- diff_pc_vs_ddd(
    cond = 1,
    seed = 2,
    pars_m = pars_m,
    pars_s = pars_s,
    age = age,
    n_0 = n_0,
    t_d = t_d,
    precision_threshold = precision_threshold
  )
  testthat::expect_true(
    test < precision_threshold || is.infinite(test)
  )

  skip("this doesn't work (yet)")

  if (packageVersion(pkg = "DDD") >= 3.8) {

    # test cond == 4
    test <- diff_pc_vs_ddd(
      cond = 4,
      seed = 2,
      pars_m = pars_m,
      pars_s = pars_s,
      age = age,
      n_0 = n_0,
      t_d = t_d,
      precision_threshold = precision_threshold
    )
    testthat::expect_true(
      test < precision_threshold || is.infinite(test)
    )
  }
})
