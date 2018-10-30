context("simulations2")

syntetic_L <- function(
  Nspecies = 30,
  age = 10,
  shifted = FALSE,
  LS
)
{
  tshift = LS[2, 1]
  Nsup <- Nspecies * 2
  L <- matrix(0, ncol = 5, nrow = Nsup)
  L[,5] <- 0
  L[,4] <- -1
  L[,3] <- c((1:Nspecies) * rep(c(1, -1), Nspecies)[1:Nspecies], rep(0, (Nsup - Nspecies)))
  L[,2] <- c(0, head(L[,3], Nsup - 1))
  L[,1] <- sort(c(age, age, runif(min = tshift, max = 10, n = Nspecies - 2), rep(0, (Nsup - Nspecies))), decreasing = TRUE)
  L[(Nspecies + 1):Nsup, 1:3] <- 0
  if (shifted == TRUE)
  {
    L <- sls_sim.use_event(
      LL = list(L),
      clade = 1,
      event = 'shift',
      t = tshift,
      LS = sls:::sls_sim.get_standard_LS()
    )$L

  }
  L
}

test_that("sls_sim.get_pars", {

  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  testthat::expect_true(
    length(pars) == length(lambdas)
  )
})

test_that("sls_sim.initialize_LL_new_clade", {

  LS <- sls:::sls_sim.get_standard_LS()
  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  Nclades <- length(lambdas)

  suppressWarnings(rm(LL))
  clade <- 1; LL <- vector("list", Nclades)
  LL <- sls:::sls_sim.initialize_LL_new_clade(LL = LL, pars = pars, clade = clade)
  pool <- sls:::sls_sim.get_pool(LL[[clade]])
  testthat::expect_true(
    length(pool) == nrow(LS)
  )
  N0s <- LS$N0
  testthat::expect_true(
    all(LL[[clade]][(1:N0s[clade]), 4] == -1)
  )

  clade <- 2
  LL <- sls_sim.initialize_LL_new_clade(pars = pars, clade = clade, LL = LL, LS = LS)
  testthat::expect_true(
    is.null(LL[[clade]])
  )
  testthat::expect_true(
    length(LL) == clade
  )
})

test_that("sls_sim.sample_deltas", {

  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf); Nclades <- length(lambdas)
  LS <- sls:::sls_sim.get_standard_LS()
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)

  suppressWarnings(rm(LL))
  clade <- 1; LL <- vector("list", Nclades)
  LL <- sls_sim.initialize_LL_new_clade(
    pars = pars,
    clade = clade,
    LL = LL,
    LS = LS
  )

  out <- sls_sim.sample_deltas(pars = pars, clade = clade, LL = LL)
  testthat::expect_true(
    out$deltaN == 1 || out$deltaN == -1
  )
  testthat::expect_true(
    out$deltaT >= 0
  )
})

test_that("sls_sim.decide_event", {

  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf); Nclades <- length(lambdas)
  LS <- sls:::sls_sim.get_standard_LS()
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  suppressWarnings(rm(LL))
  clade <- 1; LL <- vector("list", Nclades)
  LL <- sls_sim.initialize_LL_new_clade(pars = pars, clade = clade, LL = LL, LS = LS)

  t <- 0.5
  deltaN <- 1
  deltaT <- 1
  clade <- 1
  test_end <- sls_sim.decide_event(
    deltaN = deltaN,
    deltaT = deltaT,
    t = t,
    LL = LL,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_end == 'end'
  )

  t <- 5.5
  deltaN <- 1
  deltaT <- 1
  clade <- 1
  test_speciation <- sls_sim.decide_event(
    deltaN = deltaN,
    deltaT = deltaT,
    t = t,
    LL = LL,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_speciation == 'speciation'
  )

  t <- 5.5
  deltaN <- -1
  deltaT <- 1
  clade <- 1
  test_extinction <- sls_sim.decide_event(
    deltaN = deltaN,
    deltaT = deltaT,
    t = t,
    LL = LL,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_extinction == 'extinction'
  )

  # it should not shift if the shift is already saved in L
  t <- 4.5
  deltaN <- 1
  deltaT <- 1
  clade <- 1
  tshift <- 4
  LL <- sls_sim.initialize_LL_new_clade(pars = pars, clade = clade, LL = LL, LS = LS)
  LL[[clade]][3, ] <- c(6, -2, -3, tshift, 2)
  test_not_shift <- sls_sim.decide_event(
    deltaN = deltaN,
    deltaT = deltaT,
    t = t,
    LL = LL,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_not_shift != 'shift'
  )
  testthat::expect_true(
    test_not_shift == ifelse(deltaN > 0, 'speciation', 'extinction')
  )

  t <- 4.5
  deltaN <- -1
  deltaT <- 1
  clade <- 1
  LL <- sls_sim.initialize_LL_new_clade(pars = pars, clade = clade, LL = LL, LS = LS)
  LL[[clade]][3, ] <- c(6, -2, -3, -1, 0)
  test_shift <- sls_sim.decide_event(
    deltaN = deltaN,
    deltaT = deltaT,
    t = t,
    LL = LL,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_shift == 'shift'
  )
})


test_that("sls_sim.use_event", {

  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf); Nclades <- length(lambdas)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  LS <- sls:::sls_sim.get_standard_LS()
  suppressWarnings(rm(LL))
  clade <- 1; LL <- vector("list", Nclades)

  ### event 1
  event <- 'speciation'
  clade <- 1
  t <- 4.5
  L0 <- syntetic_L(Nspecies = Nspecies <- 40, LS = LS)
  LL <- sls_sim.initialize_LL_new_clade(pars = pars, clade = clade, LL = LL, LS = LS)
  LL[[clade]] <- L0
  out <- sls_sim.use_event(
    LL = LL,
    clade = clade,
    event = event,
    t = t,
    LS = sls:::sls_sim.get_standard_LS()
  )

  testthat::expect_true(
    all(L0[1:Nspecies,] == out$L[1:Nspecies,])
  )
  testthat::expect_true(
    out$L[Nspecies + 1, 1] == t
  )
  testthat::expect_true(
    abs(out$L[Nspecies + 1, 3]) == Nspecies + 1
  )
  testthat::expect_true(
    abs(out$L[Nspecies + 1, 2]) %in%
      abs(L0[1:Nspecies, 3])
  )

  ### event 2
  event <- 'extinction'
  clade <- 1
  t <- 4.5
  L0 <- syntetic_L(Nspecies = Nspecies <- 40, LS = LS)
  LL <- sls_sim.initialize_LL_new_clade(pars = pars, clade = clade, LL = LL, LS = LS)
  LL[[clade]] <- L0
  out <- sls_sim.use_event(
    LL = LL,
    clade = clade,
    event = event,
    t = t,
    LS = sls:::sls_sim.get_standard_LS()
  )

  testthat::expect_true(
    any(out$L[, 4] == t)
  )
  dead <- which(out$L[, 4] == t)
  testthat::expect_true(
    out$L[dead, 2] %in% L0[1:Nspecies, 3]
  )

  ### event 3
  event <- 'shift'
  clade <- 1
  t <- 4.5
  L0 <- syntetic_L(Nspecies = Nspecies <- 40, LS = LS)
  LL <- sls_sim.initialize_LL_new_clade(pars = pars, clade = clade, LL = LL, LS = LS)
  LL[[clade]] <- L0
  out <- sls_sim.use_event(
    LL = LL,
    clade = clade,
    event = event,
    t = t,
    LS = sls:::sls_sim.get_standard_LS()
  )

  newclade <- 2
  testthat::expect_true(
    any(out$L[1:(Nspecies + 1), 5] == newclade)
  )
  testthat::expect_true(
    out$t == LS[LS[,3]== newclade, 1]
  )

  ### event 4
  event <- 'end'
  clade <- 1
  t <- 4.5
  L0 <- syntetic_L(Nspecies = Nspecies <- 40, LS = LS)
  LL <- sls_sim.initialize_LL_new_clade(pars = pars, clade = clade, LL = LL, LS = LS)
  LL[[clade]] <- L0
  out <- sls_sim.use_event(
    LL = LL,
    clade = clade,
    event = event,
    t = t,
    LS = sls:::sls_sim.get_standard_LS()
  )

  testthat::expect_true(
    out$t == 0
  )
})

test_that("sls_sim.check_conditioning", {

  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf); Nclades <- length(lambdas)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  LS = sls:::sls_sim.get_standard_LS()
  suppressWarnings(rm(LL)); LL <- list()
  LL[[1]] <- syntetic_L(Nspecies = (Nspecies <- 5), shifted = TRUE, LS = LS); sls_sim.read_LL(LL)
  LL <- sls::sls_sim.initialize_LL_new_clade(pars = pars, clade = 2, LL = LL, LS = LS); sls_sim.read_LL(LL)
  tshift <- LS[2, 1]

  conds <- c(1, 3, 4)
  out <- rep(NA, length(conds))
  for (i in seq_along(conds))
  {
    out[i] <- sls_sim.conditioning(
      cond = conds[i],
      LL = LL,
      LS = LS
    )
  }; out

  testthat::expect_true(
    all((out <= 1) & (out >= 0) & is.numeric(out))
  )
  no_deads <- all(unique(c(
    LL[[1]][, 4],
    LL[[2]][, 4]
  )) %in% c(-1, tshift))
  all_conds_equal_one <- (all(out == 1))
  if (no_deads) {
    testthat::expect_true(
      all_conds_equal_one
    )
  }

})

test_that("sls_sim2", {

  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf); Nclades <- length(lambdas)
  LS <- sls:::sls_sim.get_standard_LS(crown_age = 5, shift_time = 2)

  maxsims <- 2
  maxtravis <- (13 * ribir::is_on_travis())
  conds <- c(3, 4)
  i <- 1
  out <- vector("list", length(conds) * (maxsims + maxtravis))
  for (s in 1:(maxsims + maxtravis)) {
    for (cond in conds) {
      if (s <= maxsims) {set.seed(s)}
      out[[i]] <- sls_sim2(
        lambdas = lambdas,
        mus = mus,
        Ks = Ks,
        cond = cond,
        LS = LS
      )

      test <- out[[i]]
      L1   <- test[[1]]
      L2   <- test[[2]]
      testthat::expect_true(
        ncol(L1) == ncol(L2),
        ncol(L1) == 5
      )
      testthat::expect_true(
        all(L1[-1, 2] %in% L1[, 3]),
        all(L2[-1, 2] %in% L2[, 3])
      )
      if (cond == 3)
      {
        testthat::expect_true(
          survM <- length(L1[L1[, 4] == -1, 3]) > 0,
          length(L2[L2[, 4] == -1, 3]) > 0
        )
        if (survM) {testthat::expect_true(sum(unique(sign(L1[L1[, 4] == -1, 3]))) == 0)}
      }
      if (cond == 4)
      {
        testthat::expect_true(
          length(L1[L1[, 4] == -1, 3]) > 0,
          length(L2[L2[, 4] == -1, 3]) > 0
        )
      }

      i <- i + 1
    }
  }

})
