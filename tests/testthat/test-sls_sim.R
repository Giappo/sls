context("sls_sim")

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

test_that("sls_sim", {

  lambdas <- c(0.2, 0.4)
  mus <- c(0.1, 0.05)
  ks <- c(Inf, Inf)
  n_clades <- length(lambdas)
  l_2 <- sls::sim_get_standard_l_2(crown_age = 5, shift_time = 2)

  time <- format(Sys.time(), "%X")
  maxsims <- 100
  maxtravis <- (
    is_on_ci() *
      (200 + 700 * (time > "23:00:00" && time < "9:30:00"))
  )
  seed_interval <- 1:(maxsims + maxtravis)
  conds <- c(3, 4)
  i <- 1
  out <- vector(
    "list",
    length(conds) * (maxsims + maxtravis)
  ); seed <- 1; cond <- conds[1]
  for (seed in seed_interval) {
    for (cond in conds) {
      if (seed <= maxsims + maxtravis / 2) {
        set.seed(seed)
      } else {
        set.seed(seed)
        lambdas <- c(x <- runif(n = 1, min = 0.2, max = 1), x / 2)
        mus <- lambdas * runif(n = 2, min = 0.2, max = 0.8)
      }
      out[[i]] <- sls_sim(
        lambdas = lambdas,
        mus = mus,
        ks = ks,
        cond = cond,
        l_2 = l_2
      )

      test <- out[[i]]
      l_0_1   <- test$l_tables[[1]]
      l_0_2   <- test$l_tables[[2]]
      testthat::expect_true(
        ncol(l_0_1) == ncol(l_0_2),
        ncol(l_0_1) == 5
      )
      testthat::expect_true(
        all(l_0_1[-1, 2] %in% l_0_1[, 3]),
        all(l_0_2[-1, 2] %in% l_0_2[, 3])
      )
      survivors_m <- l_0_1[l_0_1[, 4] == -1, 3]
      survivors_s <- l_0_2[l_0_2[, 4] == -1, 3]
      if (cond == 3) {
        testthat::expect_true(
          # at least one survivor for M and one for S
          length(survivors_m) > 0,
          length(survivors_s) > 0
        )
      }
      if (cond == 4) {
        testthat::expect_true(
          # does M survive?
          surv_m <- length(survivors_m) > 0,
          # does S survive?
          length(survivors_s) > 0
        )
        if (surv_m) {
          testthat::expect_true(
            # both left and right crown species survive
            sum(unique(sign(survivors_m))) == 0
          )
        }
      }
      i <- i + 1
    }
  }

})

test_that("sls_sim - pathological cases", {
  set.seed(1)
  l_2 <- sls::sim_get_standard_l_2(crown_age = 10, shift_time = 2)
  sim <- sls::sls_sim(
    lambdas = c(0.5399258, 0),
    mus = c(0.5400, 0),
    ks = c(Inf, Inf),
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

test_that("shift is always recorded in main clade l_0_after and
          sub clade ids always have the same sign", {

  seed_interval <- (23 - 5 * is_on_ci()):(25 + 5 * is_on_ci())
  for (seed in seed_interval) {
    set.seed(seed)
    lambdas <- c(0.3, 0.6)
    mus <- c(0.2, 0.1)
    cond <- 3
    l_2 <- sls::sim_get_standard_l_2(crown_age = 10, shift_time = 4)

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
      length(unique(sign(sim_get_pool(sim$l_tables[[2]])))) == 1
    )
  }

})

test_that("l_0 matrix size check is working", {

  lambdas <- c(0.3, 0.6)
  mus <- c(0.2, 0.1)
  cond <- 3
  l_2 <- sls::sim_get_standard_l_2(crown_age = 10, shift_time = 4)
  starting_l_size <- 10

  for (seed in 1:20) {
    set.seed(seed)
    sim <- sls::sls_sim(
      lambdas = lambdas,
      mus = mus,
      cond = cond,
      l_2 = l_2,
      l_matrix_size = starting_l_size
    )

    for (i in seq_along(nrow(l_2))) {
      testthat::expect_true(
        length(sim$brts[[i]]) <= nrow(sim$l_tables[[i]])
      )
    }
  }
})
