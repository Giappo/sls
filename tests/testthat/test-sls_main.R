context("sls_main")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("use", {
  sim_pars <- c(0.3, 0.2, 0.6, 0.1)
  cond <- 3
  loglik_functions <- sls_logliks_div()
  l_2 <- sim_get_standard_l_2(
    crown_age = 5,
    shift_time = 2
  )
  n_0 <- l_2$n_0[1]
  t_0s <- l_2$birth_time
  optim_ids <- c(TRUE, TRUE, TRUE, TRUE)
  seed_interval <- 6:(6 + 5 * is_on_ci()); seed <- seed_interval[1]
  for (seed in seed_interval) {
    # seed = 6 is critical!
    test <- sls_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      l_2 = l_2,
      start_pars = c(0.2, 0.1, 0.2, 0.1),
      loglik_functions = loglik_functions,
      optim_ids = optim_ids,
      verbose = FALSE
    )
    testthat::expect_true(
      is.data.frame(test)
    )
    testthat::expect_true(
      length(test$sim_lambda_m) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$sim_mu_m) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$sim_lambda_s) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$sim_mu_s) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$lambda_m) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$mu_m) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$lambda_s) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$mu_s) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$loglik) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$df) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$conv) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$tips_1) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$tips_2) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$seed) == length(loglik_functions)
    )
    testthat::expect_true(
      length(unique(test$seed)) == 1
    )
    testthat::expect_true(
      length(test$tips_1) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$tips_2) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$cond) == length(loglik_functions)
    )
    testthat::expect_true(
      all(
        is.numeric(test$cond)
      )
    )
    testthat::expect_true(
      all(
        test$n_0 == l_2$n_0[1]
      )
    )
    testthat::expect_true(
      all(
        test$t_0_1 == l_2$birth_time[1]
      )
    )
    testthat::expect_true(
      all(
        test$t_0_2 == l_2$birth_time[2]
      )
    )
    testthat::expect_true(
      length(test$optim_lambda_m) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$optim_mu_m) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$optim_lambda_s) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$optim_mu_s) == length(loglik_functions)
    )
    testthat::expect_true(
      all(
        c(
          test$optim_lambda_m,
          test$optim_mu_m,
          test$optim_lambda_s,
          test$optim_mu_s
        ) %in% c("TRUE", "FALSE")
      )
    )
    testthat::expect_true(
      length(test$model) == length(loglik_functions)
    )

    # test file saving
    pkg_name <- get_pkg_name() # nolint internal function
    if (.Platform$OS.type == "windows") {
      project_folder <- system.file("extdata", package = pkg_name)
    } else {
      project_folder <- getwd()
    }
    # check data folder existence
    data_folder <- file.path(project_folder, "data")
    testthat::expect_true(
      file.exists(data_folder)
    )
    # check results folder existence
    results_folder <- file.path(project_folder, "results")
    testthat::expect_true(
      file.exists(results_folder)
    )
    # check data file existence
    data_file_name <- create_data_file_name( # nolint internal function
      data_folder = data_folder,
      sim_pars = sim_pars,
      optim_ids = optim_ids,
      cond = cond,
      n_0 = n_0,
      t_0s = t_0s,
      seed = seed
    )
    testthat::expect_true(
      file.exists(data_file_name)
    )
    suppressWarnings(file.remove(data_file_name))
    # check results file existence
    results_file_name <- create_results_file_name( # nolint internal function
      results_folder = results_folder,
      sim_pars = sim_pars,
      optim_ids = optim_ids,
      cond = cond,
      n_0 = n_0,
      t_0s = t_0s,
      seed = seed
    )
    testthat::expect_true(
      file.exists(results_file_name)
    )
    # check if saved results are the right ones
    testthat::expect_equal(
      utils::read.csv(results_file_name)[, -1],
      test
    )
    suppressWarnings(file.remove(results_file_name))
  }
  # test silent mode and character entry for "loglik_functions" input
  seed <- 99
  testthat::expect_silent(
    test <- sls_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      l_2 = l_2,
      start_pars = c(0.2, 0.1, 0.2, 0.1),
      loglik_functions = "loglik_sls_p",
      optim_ids = optim_ids,
      verbose = FALSE
    )
  )
  testthat::expect_true(
    is.data.frame(test)
  )
  # check data file existence
  data_file_name <- create_data_file_name( # nolint internal function
    data_folder = data_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed
  )
  testthat::expect_true(
    file.exists(data_file_name)
  )
  suppressWarnings(file.remove(data_file_name))
  # check results file existence
  results_file_name <- create_results_file_name( # nolint internal function
    results_folder = results_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed
  )
  testthat::expect_true(
    file.exists(results_file_name)
  )
  # check if saved results are the right ones
  testthat::expect_equal(
    utils::read.csv(results_file_name)[, -1],
    test
  )
  suppressWarnings(file.remove(results_file_name))
})

test_that("it saves only once", {

  if (!is_on_ci()) {
    testthat::expect_equal(1, 1)
  } else {
    seed <- 1
    sim_pars <- c(0.27, 0.15, 0.5, 0.1)
    cond <- 3
    crown_age <- 2.5
    shift_time <- 1.8
    project_folder <- tempdir()
    results_folder <- file.path(
      project_folder,
      "results"
    )

    fn <- create_results_file_name(
      results_folder = results_folder,
      sim_pars = sim_pars,
      optim_ids = rep(TRUE, length(sim_pars)),
      cond = cond,
      n_0 = 2,
      t_0s = c(crown_age, shift_time),
      seed = seed
    ); fn

    test_p <- sls_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      start_pars = sim_pars,
      loglik_functions = loglik_sls_p,
      l_2 = sim_get_standard_l_2(
        crown_age = crown_age,
        shift_time = shift_time
      ),
      verbose = FALSE,
      project_folder = project_folder
    )
    x <- read.csv(
      file = fn
    )[, -1]

    test_p2 <- sls_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      start_pars = sim_pars,
      loglik_functions = loglik_sls_p,
      l_2 = sim_get_standard_l_2(
        crown_age = crown_age,
        shift_time = shift_time
      ),
      verbose = FALSE,
      project_folder = project_folder
    )
    y <- read.csv(
      file = fn
    )[, -1]

    testthat::expect_equal(
      x, y
    )
  }
})

test_that("it works also for a subset of parameters", {
  seed <- 10
  sim_pars <- c(0.3, 0.2, 0.6, 0.1)
  cond <- 3
  l_2 <- sim_get_standard_l_2(
    crown_age = 5,
    shift_time = 2
  )
  n_0 <- l_2$n_0[1]
  t_0s <- l_2$birth_time
  loglik_functions <- sls_logliks_div()
  optim_ids <- c(TRUE, FALSE, FALSE, FALSE)

  test <- sls_main(
    seed = seed,
    sim_pars = sim_pars,
    cond = cond,
    l_2 = l_2,
    start_pars = c(0.2, 0.1, 0.2, 0.1),
    loglik_functions = loglik_functions,
    verbose = FALSE,
    optim_ids = optim_ids
  )
  testthat::expect_true(
    is.data.frame(test)
  )
  testthat::expect_true(
    all(test$sim_mu_m == test$mu_m)
  )
  testthat::expect_true(
    all(test$sim_lambda_s == test$lambda_s)
  )
  testthat::expect_true(
    all(test$sim_mu_s == test$mu_s)
  )
  pkg_name <- get_pkg_name() # nolint internal function
  # test file saving
  if (.Platform$OS.type == "windows") {
    project_folder <- system.file("extdata", package = pkg_name)
  } else {
    project_folder <- getwd()
  }
  # check data folder existence
  data_folder <- file.path(project_folder, "data")
  testthat::expect_true(
    file.exists(data_folder)
  )
  # check results folder existence
  results_folder <- file.path(project_folder, "results")
  testthat::expect_true(
    file.exists(results_folder)
  )
  # check data file existence
  data_file_name <- create_data_file_name( # nolint internal function
    data_folder = data_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed
  )
  testthat::expect_true(
    file.exists(data_file_name)
  )
  suppressWarnings(file.remove(data_file_name))
  # check results file existence
  results_file_name <- create_results_file_name( # nolint internal function
    results_folder = results_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed
  )
  testthat::expect_true(
    file.exists(results_file_name)
  )
  suppressWarnings(file.remove(results_file_name))
})

test_that("from different likelihoods different results", {

  skip("This is long")
  if (!is_on_ci()) {skip("This test should run only on ci")}

  seed <- 5
  sim_pars <- c(0.4, 0.2, 0.6, 0.15)
  cond <- 2
  crown_age  <- 10
  shift_time <- 6
  l_2 <- sls::sim_get_standard_l_2(
    crown_age = crown_age,
    shift_time = shift_time
  )

  # test p-equation
  test <- sls_main(
    seed = seed,
    sim_pars = sim_pars,
    cond = cond,
    start_pars = sim_pars,
    loglik_functions = c(
      loglik_sls_p,
      loglik_sls_p_nodiv,
      loglik_sls_q,
      loglik_sls_q_nodiv
      ),
    l_2 = l_2,
    verbose = TRUE
  )

  testthat::expect_true(
    length(unique(unlist(
      test[, (length(sim_pars) + 1):(length(sim_pars) + 2)]
    ))) ==
      length(unlist(
        test[, (length(sim_pars) + 1):(length(sim_pars) + 2)]
      ))
  )
})

test_that("abuse", {
  seed <- 1
  sim_pars <- c(0.3, 0.2, 0.6, 0.1)
  cond <- 3
  testthat::expect_error(
    sls_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      l_2 = sim_get_standard_l_2(
        crown_age = 5,
        shift_time = 2
      ),
      start_pars = c(0.2, 0.1, 0.2, 0.1),
      loglik_functions = "nonsense",
      verbose = FALSE
    ),
    paste0(
      "This is not a likelihood function provided by ",
      get_pkg_name(),
      "!"
    )
  )
})
