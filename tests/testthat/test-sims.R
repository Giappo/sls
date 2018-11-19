context("simulations")

syntetic_data <- function(
  n_species = 30,
  shifted = FALSE,
  LS,
  clade = 1
) {
  age <- LS[1, 1]
  tshift <- LS[2, 1]
  Nsup <- n_species * 2
  L <- matrix(0, ncol = 5, nrow = Nsup)
  L[, 5] <- 0
  L[, 4] <- -1
  L[, 3] <- c(
    (1:n_species) * rep(c(1, -1), n_species)[1:n_species],
    rep(0, (Nsup - n_species))
  )
  L[, 2] <- c(0, head(L[, 3], Nsup - 1))
  L[, 1] <- sort(c(
    age,
    age,
    runif(min = tshift, max = 10, n = n_species - 2),
    rep(0, (Nsup - n_species))
  ), decreasing = TRUE)
  L[(n_species + 1):Nsup, 1:3] <- 0

  data <- list()
  data$LL[[clade]] <- L
  data$pools[[clade]] <- sls_sim.get_pool(L)
  data$Nmax <- length(
    unique(data$LL[[clade]][, 3])[unique(data$LL[[clade]][, 3]) != 0]
  )

  if (shifted == TRUE) {
    data_new <- sls_sim.use_event(
      data = data,
      clade = 1,
      event = "shift",
      t = tshift,
      LS = sls::sls_sim.get_standard_LS()
    )$data
    data <- data_new
  }
  data
}

test_that("sls_sim.get_pars", {

  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  testthat::expect_true(
    length(pars) == length(lambdas)
  )
})

test_that("sls_sim.initialize_LL_new_clade", {

  LS <- sls::sls_sim.get_standard_LS()
  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  n_clades <- length(lambdas)

  suppressWarnings(rm(data))
  data <- sls::sls_sim.initialize_LL_new_clade(
    data = data,
    pars = pars,
    clade = 0,
    LS = LS
  )

  clade <- 1
  data <- sls::sls_sim.initialize_LL_new_clade(
    data = data,
    pars = pars,
    clade = clade,
    LS = LS
  )
  testthat::expect_true(
    length(data$pools[[clade]]) == nrow(LS)
  )
  n_0s <- LS$n_0
  testthat::expect_true(
    all(data$LL[[clade]][(1:n_0s[clade]), 4] == -1)
  )

  clade <- 2
  data <- sls_sim.initialize_LL_new_clade(
    data = data,
    pars = pars,
    clade = clade,
    LS = LS
  )
  testthat::expect_true(
    is.null(data$LL[[clade]])
  )
  testthat::expect_true(
    length(data$LL) == clade
  )
})

test_that("sls_sim.sample_deltas", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  Ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  LS <- sls::sls_sim.get_standard_LS()
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)

  suppressWarnings(rm(data))
  clade <- 1; data <- sls_sim.initialize_LL_new_clade(pars = pars, clade = 0)
  data <- sls_sim.initialize_LL_new_clade(
    data = data,
    pars = pars,
    clade = clade,
    LS = LS
  )

  out <- sls_sim.sample_deltas(pars = pars, clade = clade, data = data)
  testthat::expect_true(
    out$delta_n == 1 || out$delta_n == -1
  )
  testthat::expect_true(
    out$delta_t >= 0
  )
})

test_that("sls_sim.decide_event", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  Ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  LS <- sls::sls_sim.get_standard_LS()
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  suppressWarnings(rm(LL))
  clade <- 1; data <- sls_sim.initialize_LL_new_clade(pars = pars, clade = 0)
  data <- sls_sim.initialize_LL_new_clade(
    pars = pars,
    clade = clade,
    data = data,
    LS = LS
  )

  t <- 0.5
  delta_n <- 1
  delta_t <- 1
  clade <- 1
  test_end <- sls_sim.decide_event(
    delta_n = delta_n,
    delta_t = delta_t,
    t = t,
    data = data,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_end == "end"
  )

  t <- 5.5
  delta_n <- 1
  delta_t <- 1
  clade <- 1
  test_speciation <- sls_sim.decide_event(
    delta_n = delta_n,
    delta_t = delta_t,
    t = t,
    data = data,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_speciation == "speciation"
  )

  t <- 5.5
  delta_n <- -1
  delta_t <- 1
  clade <- 1
  test_extinction <- sls_sim.decide_event(
    delta_n = delta_n,
    delta_t = delta_t,
    t = t,
    data = data,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_extinction == "extinction"
  )

  # it should not shift if the shift is already saved in L
  t <- 4.5
  delta_n <- 1
  delta_t <- 1
  clade <- 1
  tshift <- 4
  data <- sls_sim.initialize_LL_new_clade(
    pars = pars,
    clade = clade,
    data = data,
    LS = LS
  )
  data$LL[[clade]][3, ] <- c(6, -2, -3, tshift, 2)
  data$Nmax <- length(
    unique(data$LL[[clade]][, 3])[unique(data$LL[[clade]][, 3]) != 0]
  )
  test_not_shift <- sls_sim.decide_event(
    delta_n = delta_n,
    delta_t = delta_t,
    t = t,
    data = data,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_not_shift != "shift"
  )
  testthat::expect_true(
    test_not_shift == ifelse(delta_n > 0, "speciation", "extinction")
  )

  t <- 4.5
  delta_n <- -1
  delta_t <- 1
  clade <- 1
  data <- sls_sim.initialize_LL_new_clade(
    pars = pars,
    clade = clade,
    data = data,
    LS = LS
  )
  data$LL[[clade]][3, ] <- c(6, -2, -3, -1, 0)
  test_shift <- sls_sim.decide_event(
    delta_n = delta_n,
    delta_t = delta_t,
    t = t,
    data = data,
    clade = clade,
    LS = LS
  )

  testthat::expect_true(
    test_shift == "shift"
  )
})


test_that("sls_sim.use_event", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  Ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  LS <- sls::sls_sim.get_standard_LS()
  suppressWarnings(rm(data))
  data <- sls_sim.initialize_LL_new_clade(pars = pars, clade = 0)

  ### event 1
  event <- "speciation"
  clade <- 1
  t <- 4.5
  data <- syntetic_data(n_species = n_species <- 40, LS = LS)
  L0 <- data$LL[[clade]]
  out <- sls_sim.use_event(
    data = data,
    clade = clade,
    event = event,
    t = t,
    LS = sls::sls_sim.get_standard_LS()
  ); L <- out$data$LL[[clade]]

  testthat::expect_true(
    all(
      L0[1:n_species, ] ==
      L[1:n_species, ]
    )
  )
  testthat::expect_true(
    L[n_species + 1, 1] == t
  )
  testthat::expect_true(
    abs(L[n_species + 1, 3]) == n_species + 1
  )
  testthat::expect_true(
    abs(L[n_species + 1, 2]) %in%
      abs(L0[1:n_species, 3])
  )

  ### event 2
  event <- "extinction"
  clade <- 1
  t <- 4.5
  data <- syntetic_data(n_species = n_species <- 40, LS = LS)
  L0 <- data$LL[[clade]]
  out <- sls_sim.use_event(
    data = data,
    clade = clade,
    event = event,
    t = t,
    LS = sls::sls_sim.get_standard_LS()
  ); L <- out$data$LL[[clade]]

  testthat::expect_true(
    any(L[, 4] == t)
  )
  dead <- which(L[, 4] == t)
  testthat::expect_true(
    L[dead, 2] %in% c(0, L0[1:n_species, 3])
  )

  ### event 3
  event <- "shift"
  clade <- 1
  t <- 4.5
  data <- syntetic_data(n_species = n_species <- 40, LS = LS)
  L0 <- data$LL[[clade]]
  out <- sls_sim.use_event(
    data = data,
    clade = clade,
    event = event,
    t = t,
    LS = sls::sls_sim.get_standard_LS()
  ); L <- out$data$LL[[clade]]

  newclade <- LS$clade_id[2]
  testthat::expect_true(
    any(L[1:(n_species + 1), 5] == newclade)
  )
  testthat::expect_true(
    out$t == LS[LS[, 3] == newclade, 1]
  )

  ### event 4
  event <- "end"
  clade <- 1
  t <- 4.5
  data <- syntetic_data(n_species = n_species <- 40, LS = LS)
  L0 <- data$LL[[clade]]
  out <- sls_sim.use_event(
    data = data,
    clade = clade,
    event = event,
    t = t,
    LS = sls::sls_sim.get_standard_LS()
  ); L <- out$data$LL[[clade]]

  testthat::expect_true(
    out$t == 0
  )
})

test_that("sls_sim.check_conditioning", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  Ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  LS <- sls::sls_sim.get_standard_LS()
  suppressWarnings(rm(data));
  data <- syntetic_data(n_species = (n_species <- 5), shifted = TRUE, LS = LS);
  data <- sls::sls_sim.initialize_LL_new_clade(
    pars = pars,
    clade = 2,
    data = data,
    LS = LS
  ); sls_sim.read_LL(data$LL)
  tshift <- LS[2, 1]

  conds <- c(1, 3, 4)
  out <- rep(NA, length(conds))
  for (i in seq_along(conds)) {
    out[i] <- sls_sim.conditioning(
      cond = conds[i],
      data = data,
      LS = LS
    )
  }

  testthat::expect_true(
    all(
      (out <= 1) & (out >= 0) & is.numeric(out)
    )
  )
  no_deads <- all(unique(c(
    data$LL[[1]][, 4],
    data$LL[[2]][, 4]
  )) %in% c(-1, tshift))
  all_conds_equal_one <- (all(out == 1))
  if (no_deads) {
    testthat::expect_true(
      all_conds_equal_one
    )
  }

})

test_that("sls_sim", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  Ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  LS <- sls::sls_sim.get_standard_LS(crown_age = 5, shift_time = 2)

  maxsims <- 2
  maxtravis <- (13 * ribir::is_on_travis())
  conds <- c(3, 4)
  i <- 1
  out <- vector(
    "list",
    length(conds) * (maxsims + maxtravis)
  ); s <- 1; cond <- conds[1]
  for (s in 1:(maxsims + maxtravis)) {
    for (cond in conds) {
      if (s <= maxsims) {
        set.seed(s)
      }
      out[[i]] <- sls_sim(
        lambdas = lambdas,
        mus = mus,
        Ks = Ks,
        cond = cond,
        LS = LS
      )

      test <- out[[i]]
      L1   <- test$l_tables[[1]]
      L2   <- test$l_tables[[2]]
      testthat::expect_true(
        ncol(L1) == ncol(L2),
        ncol(L1) == 5
      )
      testthat::expect_true(
        all(L1[-1, 2] %in% L1[, 3]),
        all(L2[-1, 2] %in% L2[, 3])
      )
      if (cond == 3) {
        testthat::expect_true(
          surv_m <- length(L1[L1[, 4] == -1, 3]) > 0,
          length(L2[L2[, 4] == -1, 3]) > 0
        )
        if (surv_m) {
          testthat::expect_true(sum(unique(sign(L1[L1[, 4] == -1, 3]))) == 0)
        }
      }
      if (cond == 4) {
        testthat::expect_true(
          length(L1[L1[, 4] == -1, 3]) > 0,
          length(L2[L2[, 4] == -1, 3]) > 0
        )
      }

      i <- i + 1
    }
  }

})

test_that("sls_sim - pathological cases", {
  set.seed(1)
  sls::sls_sim(
    lambdas = c(0.5399258, 0),
    mus = c(0.5400, 0),
    Ks = c(Inf, Inf),
    cond = 3,
    LS = sls::sls_sim.get_standard_LS(crown_age = 10, shift_time = 2)
  )
})

test_that("shift is always recorded in main clade L and
          sub clade ids always have the same sign", {

  for (s in 23:25) {
    set.seed(s)
    print(s)
    lambdas <- c(0.3, 0.6)
    mus <- c(0.2, 0.1)
    cond <- 3
    LS <- sls::sls_sim.get_standard_LS(crown_age = 10, shift_time = 4)

    sim <- sls::sls_sim(
      lambdas = lambdas,
      mus = mus,
      cond = cond,
      LS = LS
    )
    testthat::expect_true(
      # if subclade is not created there should be
      # no sign of it in the main clade table
      !is.null(sim$l_tables[[2]]) == any(sim$l_tables[[1]][, 5] == 2)
    )
    testthat::expect_true(
      # in the subclade all the ids have the same sign
      length(unique(sign(sls_sim.get_pool(sim$l_tables[[2]])))) == 1
    )
  }

})

test_that("L matrix size check is working", {

  lambdas <- c(0.3, 0.6)
  mus <- c(0.2, 0.1)
  cond <- 3
  LS <- sls::sls_sim.get_standard_LS(crown_age = 10, shift_time = 4)
  starting_l_size <- 10

  for (s in 1:30) {
    set.seed(s)
    print(s)
    sim <- sls::sls_sim(
      lambdas = lambdas,
      mus = mus,
      cond = cond,
      LS = LS,
      l_matrix_size = starting_l_size
    ); sim

    for (i in seq_along(nrow(LS))) {
      testthat::expect_true(
        length(sim$brts[[i]]) <= nrow(sim$l_tables[[i]])
      )
    }
  }

})
