context("sls_sim_utils")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

syntetic_data <- function(
  n_species = 30,
  shifted = FALSE,
  l_2 = sls::sim_get_standard_l_2(),
  clade = 1
) {
  age <- l_2[1, 1]
  tshift <- l_2[2, 1]
  n_sup <- n_species * 2
  l_0 <- matrix(0, ncol = 5, nrow = n_sup)
  l_0[, 5] <- 0
  l_0[, 4] <- -1
  l_0[, 3] <- c(
    (1:n_species) * rep(c(1, -1), n_species)[1:n_species],
    rep(0, (n_sup - n_species))
  )
  l_0[, 2] <- c(0, head(l_0[, 3], n_sup - 1))
  l_0[, 1] <- sort(c(
    age,
    age,
    runif(min = tshift, max = 10, n = n_species - 2),
    rep(0, (n_sup - n_species))
  ), decreasing = TRUE)
  l_0[(n_species + 1):n_sup, 1:3] <- 0

  data <- list()
  data$l_1[[clade]] <- l_0
  data$pools[[clade]] <- sim_get_pool(l_0) # nolint internal function
  data$n_max <- length(
    unique(data$l_1[[clade]][, 3])[unique(data$l_1[[clade]][, 3]) != 0]
  )

  if (shifted == TRUE) {
    deltas <- list(
      delta_n = 1,
      delta_t = 0.1
    )
    data_new <- sim_use_event(
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

test_that("sim_get_pars", {

  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); ks <- c(Inf, Inf)
  pars <- sim_get_pars(lambdas = lambdas, mus = mus, ks = ks)
  testthat::expect_true(
    length(pars) == length(lambdas)
  )
})

test_that("sim_initialize_data_new_clade", {

  l_2 <- sls::sim_get_standard_l_2()
  lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); ks <- c(Inf, Inf)
  pars <- sim_get_pars(lambdas = lambdas, mus = mus, ks = ks)
  n_clades <- length(lambdas)

  suppressWarnings(rm(data))
  data <- sls::sim_initialize_data_new_clade(
    data = data,
    pars = pars,
    clade = 0,
    l_2 = l_2
  )

  clade <- 1
  data <- sls::sim_initialize_data_new_clade(
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
  data <- sim_initialize_data_new_clade(
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

  # it works even if you specify the wrong matrix size
  clade <- 1
  data <- sls::sim_initialize_data_new_clade(
    data = data,
    pars = pars,
    clade = 0,
    l_2 = l_2,
    l_matrix_size = 2
  )
  testthat::expect_true(
    length(data$l_1) > 0
  )
  data <- sim_initialize_data_new_clade(
    data = data,
    pars = pars,
    clade = clade,
    l_2 = l_2,
    l_matrix_size = 2
  )
  testthat::expect_true(
    length(data$l_1) > 0
  )
})

test_that("sim_sample_deltas", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  l_2 <- sls::sim_get_standard_l_2()
  pars <- sim_get_pars(lambdas = lambdas, mus = mus, ks = ks)

  suppressWarnings(rm(data))
  clade <- 1; data <- sim_initialize_data_new_clade(pars = pars, clade = 0)
  data <- sim_initialize_data_new_clade(
    data = data,
    pars = pars,
    clade = clade,
    l_2 = l_2
  )

  out <- sim_sample_deltas(pars = pars, clade = clade, data = data)
  testthat::expect_true(
    out$delta_n == 1 || out$delta_n == -1
  )
  testthat::expect_true(
    out$delta_t >= 0
  )
})

test_that("sim_decide_event", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  l_2 <- sls::sim_get_standard_l_2()
  pars <- sim_get_pars(lambdas = lambdas, mus = mus, ks = ks)
  suppressWarnings(rm(l_1))
  clade <- 1; data <- sim_initialize_data_new_clade(pars = pars, clade = 0)
  data <- sim_initialize_data_new_clade(
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
  test_end <- sim_decide_event(
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
  test_speciation <- sim_decide_event(
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
  test_extinction <- sim_decide_event(
    deltas = deltas,
    data = data,
    clade = clade,
    l_2 = l_2
  )

  testthat::expect_true(
    test_extinction == "extinction"
  )

  # it should not shift if the shift is already saved in l_0
  clade <- 1
  delta_n <- -1
  delta_t <- 1
  deltas <- list(
    delta_n = unname(delta_n),
    delta_t = unname(delta_t)
  )
  tshift <- 4
  data <- sim_initialize_data_new_clade(
    pars = pars,
    clade = clade,
    data = data,
    l_2 = l_2
  )
  data$l_1[[clade]][3, ] <- c(6, -2, -3, tshift, 2) # register the shift
  data$n_max[[clade]] <- length(
    unique(data$l_1[[clade]][, 3])[unique(data$l_1[[clade]][, 3]) != 0]
  )
  data$t[[clade]] <- 4.5
  test_not_shift <- sim_decide_event(
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
  data <- sim_initialize_data_new_clade(
    pars = pars,
    clade = clade,
    data = data,
    l_2 = l_2
  )
  data$l_1[[clade]][3, ] <- c(6, -2, -3, -1, 0)
  data$t[[clade]] <- 4.5
  test_shift <- sim_decide_event(
    deltas = deltas,
    data = data,
    clade = clade,
    l_2 = l_2
  )

  testthat::expect_true(
    test_shift == "shift"
  )
})


test_that("sim_use_event", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  pars <- sim_get_pars(lambdas = lambdas, mus = mus, ks = ks)
  l_2 <- sls::sim_get_standard_l_2()
  suppressWarnings(rm(data))
  data <- sim_initialize_data_new_clade(pars = pars, clade = 0)

  ### event 1
  event <- "speciation"
  n_species <- 40
  data <- syntetic_data(l_2 = l_2, n_species = n_species)
  clade <- 1
  deltas <- list(
    delta_n = 1,
    delta_t = 1
  )
  l_0_before <- data$l_1[[clade]]
  out <- sim_use_event(
    data = data,
    clade = clade,
    event = event,
    deltas = deltas,
    l_2 = l_2
  ); l_0_after <- out$l_1[[clade]]

  testthat::expect_true(
    all(
      l_0_before[1:n_species, ] ==
      l_0_after[1:n_species, ]
    )
  )
  testthat::expect_true(
    l_0_after[n_species + 1, 1] == out$t[[clade]]
  )
  testthat::expect_true(
    abs(l_0_after[n_species + 1, 3]) == n_species + 1
  )
  testthat::expect_true(
    abs(l_0_after[n_species + 1, 2]) %in%
      abs(l_0_before[1:n_species, 3])
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
  l_0_before <- data$l_1[[clade]]
  out <- sim_use_event(
    data = data,
    clade = clade,
    event = event,
    deltas = deltas,
    l_2 = l_2
  ); l_0_after <- out$l_1[[clade]]

  testthat::expect_true(
    any(l_0_after[, 4] == out$t[[clade]])
  )
  dead <- which(l_0_after[, 4] == out$t[[clade]])
  testthat::expect_true(
    l_0_after[dead, 2] %in% c(0, l_0_before[1:n_species, 3])
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
  l_0_before <- data$l_1[[clade]]
  out <- sim_use_event(
    data = data,
    clade = clade,
    event = event,
    deltas = deltas,
    l_2 = l_2
  ); l_0_after <- out$l_1[[clade]]

  newclade <- l_2$clade_id[2]
  testthat::expect_true(
    any(l_0_after[1:(n_species + 1), 5] == newclade)
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
  l_0_before <- data$l_1[[clade]]
  out <- sim_use_event(
    data = data,
    clade = clade,
    event = event,
    deltas = deltas,
    l_2 = l_2
  ); l_0_after <- out$l_1[[clade]]

  testthat::expect_true(
    out$t == 0
  )
})

test_that("sim_check_conditioning", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  pars <- sim_get_pars(lambdas = lambdas, mus = mus, ks = ks)
  l_2 <- sls::sim_get_standard_l_2()
  suppressWarnings(rm(data));
  n_species <- 5
  data <- syntetic_data(n_species = n_species, shifted = TRUE, l_2 = l_2);
  data <- sls::sim_initialize_data_new_clade(
    pars = pars,
    clade = 2,
    data = data,
    l_2 = l_2
  ); sim_read_l_1(data$l_1)
  tshift <- l_2[2, 1]

  conds <- c(1, 3, 4)
  out <- rep(NA, length(conds))
  for (i in seq_along(conds)) {
    out[i] <- sim_conditioning(
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

# Different function that should do the same thing. Used to compare.
sim_get_brts2 <- function(
  data,
  l_2
) {
  if (!all(data$t == 0)) {
    stop("times in all clades need to reach 0 before you call this function")
  }
  brts <- vector("list", length(data$l_1))

  for (clade in seq_along(data$l_1)) {
    n_0 <- l_2$n_0[clade]
    done <- 0
    if (is.null(data$l_1[[clade]]) && done == 0) {
      brts[clade] <- NULL
      done <- 1
    }
    if (done == 0) {
      l_0 <- sim_cut_l_matrix(
        unname(data$l_1[[clade]])
      )
    }
    l_0[l_0[, 5] > 0, 4] <- -1
    if (nrow(l_0) == 1 && l_0[, 4] == -1) {
      brts[[clade]] <- l_0[1, 1]
    } else {
      brts[[clade]] <- DDD::L2brts(L = l_0)
      if (n_0 == 1) {
        brts[[clade]] <- c(l_0[1, 1], brts[[clade]])
      }
    }
    brts[[clade]] <- unname(brts[[clade]])
  }
  brts <- unname(brts)
  brts
}

# Test whether the function yields the right amount of tips and
# whether these branching times are actually included in the data
test_brts_fun <- function(
  l_2,
  data,
  brts_fun = sim_get_brts
) {
  max_clade <- 2
  out <- rep(0, max_clade)
  brts <- brts_fun(
    data = data,
    l_2 = l_2
  )
  for (clade in 1:max_clade) {
    n_0 <- l_2$n_0[clade]
    br_ts <- brts[[clade]]
    tips <- (n_0 - 1) + length(br_ts); tips
    exp_tips <- (n_0 - 1) + sum(data$l_1[[clade]][, 4] == -1)
    out[clade] <- (tips == exp_tips) *
      all(
        signif(br_ts, digits = 5) %in%
          signif(data$l_1[[clade]][, 1], digits = 5)
      )
  }
  prod(out)
}

# Simulate the process until the step when the conversion from data to brts
# is needed
sim_data <- function(
  lambdas,
  mus,
  ks = c(Inf, Inf),
  l_2 = sim_get_standard_l_2(
    crown_age = 5,
    shift_time = 2
  ),
  cond,
  l_matrix_size = 1e4
) {
 # create the parameters
  pars <- sim_get_pars(
    lambdas = lambdas,
    mus = mus,
    ks = ks
  )
  good_sim <- 0
  while (!good_sim) {

    # initialize data
    data <- sim_initialize_data_new_clade(clade = 0, l_2 = l_2); clade <- 1; # nolint internal function
    for (clade in l_2$clade_id) {

      # initialize data for the clade
      data <- sim_initialize_data_new_clade(
        data = data,
        clade = clade,
        pars = pars,
        l_2 = l_2,
        l_matrix_size = l_matrix_size
      )
      while (data$t[[clade]] > 0) {

        # sample delta_n and delta_t
        deltas <- sim_sample_deltas(
          data = data,
          clade = clade,
          pars = pars
        ); deltas

        # decide the event
        event <- sim_decide_event(
          data = data,
          clade = clade,
          l_2 = l_2,
          deltas = deltas
        ); event

        # modify data accordingly
        output <- sim_use_event(
          data = data,
          clade = clade,
          l_2 = l_2,
          event = event,
          deltas = deltas
        ); output
        data <- output
      }
    }

    # is the simulation in agreement with the conditioning?
    good_sim <- sim_conditioning(
      data = data,
      l_2 = l_2,
      cond = cond
    ); good_sim
  }
  data
}

test_that("sim_get_brts", {
  lambdas <- c(0.3, 0.6)
  mus <- c(0.2, 0.1)
  l_2 <- sim_get_standard_l_2(
    crown_age = 5,
    shift_time = 2
  )
  seed_interval <- 1:20
  brts_funs <- c(
    sim_get_brts,
    sim_get_brts2
  )
  cond <- sls_conds()[1]
  brts_list <- outs <- vector("list", length(brts_funs))
  for (seed in seed_interval) {
    set.seed(seed)
    cond <- sls_conds()[1] * (cond == sls_conds()[2]) +
      sls_conds()[2] * (cond == sls_conds()[1])
    data <- sim_data(
      lambdas = lambdas,
      mus = mus,
      l_2 = l_2,
      cond = cond
    )
    for (i in seq_along(brts_funs)) {
      outs[[i]][seed] <- test_brts_fun(
        brts_fun = brts_funs[[i]],
        l_2 = l_2,
        data = data
      )
      testthat::expect_true(
        outs[[i]][seed] == 1
      )
      brts_list[[i]][[seed]] <- brts_funs[[i]](
        l_2 = l_2,
        data = data
      )
    }
  }
  for (seed in seed_interval) {
    testthat::expect_true(
      all.equal(
        brts_list[[1]][[seed]], brts_list[[2]][[seed]]
      )
    )
  }
})
