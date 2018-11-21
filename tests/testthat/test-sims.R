context("simulations")

syntetic_data <- function(
  n_species = 30,
  shifted = FALSE,
  l_2 = sls::sls_sim.get_standard_l_2(),
  clade = 1
) {
  age <- l_2[1, 1]
  tshift <- l_2[2, 1]
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
  data$l_1[[clade]] <- L
  data$pools[[clade]] <- sls_sim.get_pool(L)
  data$Nmax <- length(
    unique(data$l_1[[clade]][, 3])[unique(data$l_1[[clade]][, 3]) != 0]
  )

  if (shifted == TRUE) {
    deltas <- list(
      delta_n = 1,
      delta_t = 0.1
    )
    data_new <- sls_sim.use_event(
      data = data,
      clade = 1,
      event = "shift",
      deltas = deltas,
      l_2 = l_2
    )
    data <- data_new
  }

  event_times <- c(
    data$l_1[[clade]][data$l_1[[clade]][, 1] != 0, 1],
    data$l_1[[clade]][data$l_1[[clade]][, 4] != -1, 4]
  )
  data$t[[clade]] <- min(
    event_times
  )
  data
}

test_that("sls_sim.get_pars", {

  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  testthat::expect_true(
    length(pars) == length(lambdas)
  )
})

test_that("sls_sim.initialize_data_new_clade", {

  l_2 <- sls::sls_sim.get_standard_l_2()
  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf)
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  n_clades <- length(lambdas)

  suppressWarnings(rm(data))
  data <- sls::sls_sim.initialize_data_new_clade(
    data = data,
    pars = pars,
    clade = 0,
    l_2 = l_2
  )

  clade <- 1
  data <- sls::sls_sim.initialize_data_new_clade(
    data = data,
    pars = pars,
    clade = clade,
    l_2 = l_2
  )
  testthat::expect_true(
    length(data$pools[[clade]]) == nrow(l_2)
  )
  n_0s <- l_2$n_0
  testthat::expect_true(
    all(data$l_1[[clade]][(1:n_0s[clade]), 4] == -1)
  )

  clade <- 2
  data <- sls_sim.initialize_data_new_clade(
    data = data,
    pars = pars,
    clade = clade,
    l_2 = l_2
  )
  testthat::expect_true(
    is.null(data$l_1[[clade]])
  )
  testthat::expect_true(
    length(data$l_1) == clade
  )
})

test_that("sls_sim.sample_deltas", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  Ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  l_2 <- sls::sls_sim.get_standard_l_2()
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)

  suppressWarnings(rm(data))
  clade <- 1; data <- sls_sim.initialize_data_new_clade(pars = pars, clade = 0)
  data <- sls_sim.initialize_data_new_clade(
    data = data,
    pars = pars,
    clade = clade,
    l_2 = l_2
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
  l_2 <- sls::sls_sim.get_standard_l_2()
  pars <- sls_sim.get_pars(lambdas = lambdas, mus = mus, Ks = Ks)
  suppressWarnings(rm(l_1))
  clade <- 1; data <- sls_sim.initialize_data_new_clade(pars = pars, clade = 0)
  data <- sls_sim.initialize_data_new_clade(
    pars = pars,
    clade = clade,
    data = data,
    l_2 = l_2
  )

  ###
  clade <- 1
  delta_n <- 1
  delta_t <- 1
  deltas <- list(
    delta_n = unname(delta_n),
    delta_t = unname(delta_t)
  )
  data$t[[clade]] <- 0.5
  test_end <- sls_sim.decide_event(
    deltas = deltas,
    data = data,
    clade = clade,
    l_2 = l_2
  )

  testthat::expect_true(
    test_end == "end"
  )

  ###
  clade <- 1
  delta_n <- 1
  delta_t <- 1
  deltas <- list(
    delta_n = unname(delta_n),
    delta_t = unname(delta_t)
  )
  data$t[[clade]] <- 5.5
  test_speciation <- sls_sim.decide_event(
    deltas = deltas,
    data = data,
    clade = clade,
    l_2 = l_2
  )

  testthat::expect_true(
    test_speciation == "speciation"
  )

  ###
  clade <- 1
  delta_n <- -1
  delta_t <- 1
  deltas <- list(
    delta_n = unname(delta_n),
    delta_t = unname(delta_t)
  )
  data$t[[clade]] <- 5.5
  test_extinction <- sls_sim.decide_event(
    deltas = deltas,
    data = data,
    clade = clade,
    l_2 = l_2
  )

  testthat::expect_true(
    test_extinction == "extinction"
  )

  # it should not shift if the shift is already saved in L
  clade <- 1
  delta_n <- -1
  delta_t <- 1
  deltas <- list(
    delta_n = unname(delta_n),
    delta_t = unname(delta_t)
  )
  tshift <- 4
  data <- sls_sim.initialize_data_new_clade(
    pars = pars,
    clade = clade,
    data = data,
    l_2 = l_2
  )
  data$l_1[[clade]][3, ] <- c(6, -2, -3, tshift, 2) # register the shift
  data$Nmax[[clade]] <- length(
    unique(data$l_1[[clade]][, 3])[unique(data$l_1[[clade]][, 3]) != 0]
  )
  data$t[[clade]] <- 4.5
  test_not_shift <- sls_sim.decide_event(
    deltas = deltas,
    data = data,
    clade = clade,
    l_2 = l_2
  )

  testthat::expect_true(
    test_not_shift != "shift"
  )
  testthat::expect_true(
    test_not_shift == ifelse(delta_n > 0, "speciation", "extinction")
  )

  ###
  clade <- 1
  delta_n <- -1
  delta_t <- 1
  deltas <- list(
    delta_n = unname(delta_n),
    delta_t = unname(delta_t)
  )
  data <- sls_sim.initialize_data_new_clade(
    pars = pars,
    clade = clade,
    data = data,
    l_2 = l_2
  )
  data$l_1[[clade]][3, ] <- c(6, -2, -3, -1, 0)
  data$t[[clade]] <- 4.5
  test_shift <- sls_sim.decide_event(
    deltas = deltas,
    data = data,
    clade = clade,
    l_2 = l_2
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
  l_2 <- sls::sls_sim.get_standard_l_2()
  suppressWarnings(rm(data))
  data <- sls_sim.initialize_data_new_clade(pars = pars, clade = 0)

  ### event 1
  event <- "speciation"
  n_species <- 40
  data <- syntetic_data(l_2 = l_2, n_species = n_species)
  clade <- 1
  deltas <- list(
    delta_n = 1,
    delta_t = 1
  )
  L0 <- data$l_1[[clade]]
  out <- sls_sim.use_event(
    data = data,
    clade = clade,
    event = event,
    deltas = deltas,
    l_2 = l_2
  ); L <- out$l_1[[clade]]

  testthat::expect_true(
    all(
      L0[1:n_species, ] ==
      L[1:n_species, ]
    )
  )
  testthat::expect_true(
    L[n_species + 1, 1] == out$t[[clade]]
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
  n_species <- 40
  data <- syntetic_data(l_2 = l_2, n_species = n_species)
  clade <- 1
  deltas <- list(
    delta_n = -1,
    delta_t = 1
  )
  L0 <- data$l_1[[clade]]
  out <- sls_sim.use_event(
    data = data,
    clade = clade,
    event = event,
    deltas = deltas,
    l_2 = l_2
  ); L <- out$l_1[[clade]]

  testthat::expect_true(
    any(L[, 4] == out$t[[clade]])
  )
  dead <- which(L[, 4] == out$t[[clade]])
  testthat::expect_true(
    L[dead, 2] %in% c(0, L0[1:n_species, 3])
  )

  ### event 3
  event <- "shift"
  n_species <- 40
  data <- syntetic_data(l_2 = l_2, n_species = n_species)
  clade <- 1
  deltas <- list(
    delta_n = -1,
    delta_t = 1
  )
  L0 <- data$l_1[[clade]]
  out <- sls_sim.use_event(
    data = data,
    clade = clade,
    event = event,
    deltas = deltas,
    l_2 = l_2
  ); L <- out$l_1[[clade]]

  newclade <- l_2$clade_id[2]
  testthat::expect_true(
    any(L[1:(n_species + 1), 5] == newclade)
  )
  testthat::expect_true(
    out$t == l_2[l_2[, 3] == newclade, 1]
  )

  ### event 4
  event <- "end"
  n_species <- 40
  data <- syntetic_data(l_2 = l_2, n_species = n_species)
  clade <- 1
  deltas <- list(
    delta_n = -1,
    delta_t = 1
  )
  L0 <- data$l_1[[clade]]
  out <- sls_sim.use_event(
    data = data,
    clade = clade,
    event = event,
    deltas = deltas,
    l_2 = l_2
  ); L <- out$l_1[[clade]]

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
  l_2 <- sls::sls_sim.get_standard_l_2()
  suppressWarnings(rm(data));
  n_species <- 5
  data <- syntetic_data(n_species = n_species, shifted = TRUE, l_2 = l_2);
  data <- sls::sls_sim.initialize_data_new_clade(
    pars = pars,
    clade = 2,
    data = data,
    l_2 = l_2
  ); sls_sim.read_l_1(data$l_1)
  tshift <- l_2[2, 1]

  conds <- c(1, 3, 4)
  out <- rep(NA, length(conds))
  for (i in seq_along(conds)) {
    out[i] <- sls_sim.conditioning(
      cond = conds[i],
      data = data,
      l_2 = l_2
    )
  }

  testthat::expect_true(
    all(
      (out <= 1) & (out >= 0) & is.numeric(out)
    )
  )
  no_deads <- all(unique(c(
    data$l_1[[1]][, 4],
    data$l_1[[2]][, 4]
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
  l_2 <- sls::sls_sim.get_standard_l_2(crown_age = 5, shift_time = 2)

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
        l_2 = l_2
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
  l_2 <- sls::sls_sim.get_standard_l_2(crown_age = 10, shift_time = 2)
  sim <- sls::sls_sim(
    lambdas = c(0.5399258, 0),
    mus = c(0.5400, 0),
    Ks = c(Inf, Inf),
    cond = 3,
    l_2 = l_2
  )
  n_clades <- nrow(l_2)
  testthat::expect_true(
    length(sim$l_tables) == n_clades
  )
  testthat::expect_true(
    length(sim$brts) == n_clades
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
    l_2 <- sls::sls_sim.get_standard_l_2(crown_age = 10, shift_time = 4)

    sim <- sls::sls_sim(
      lambdas = lambdas,
      mus = mus,
      cond = cond,
      l_2 = l_2
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
  l_2 <- sls::sls_sim.get_standard_l_2(crown_age = 10, shift_time = 4)
  starting_l_size <- 10

  for (s in 1:30) {
    set.seed(s)
    print(s)
    sim <- sls::sls_sim(
      lambdas = lambdas,
      mus = mus,
      cond = cond,
      l_2 = l_2,
      l_matrix_size = starting_l_size
    ); sim

    for (i in seq_along(nrow(l_2))) {
      testthat::expect_true(
        length(sim$brts[[i]]) <= nrow(sim$l_tables[[i]])
      )
    }
  }

})
