context("sim_get_brts")

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
    data <- sim_initialize_data_new_clade(clade = 0, l_2 = l_2); clade <- 1;
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
