context("sls_main")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("use", {
  seed_interval <- 6:(6 + 5 * is_on_ci()) # 6 is critical
  for (seed in seed_interval) {
    sim_pars <- c(0.3, 0.2, 0.6, 0.1)
    cond <- 3
    models <- sls_logliks_div()
    l_2 <- sim_get_standard_l_2(
      crown_age = 5,
      shift_time = 2
    )
    test <- sls_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      l_2 = l_2,
      start_pars = c(0.2, 0.1, 0.2, 0.1),
      models = models,
      verbose = FALSE
    )
    testthat::expect_true(
      is.data.frame(test)
    )
    testthat::expect_true(
      length(test$sim_lambda_m) == length(models)
    )
    testthat::expect_true(
      length(test$sim_mu_m) == length(models)
    )
    testthat::expect_true(
      length(test$sim_lambda_s) == length(models)
    )
    testthat::expect_true(
      length(test$sim_mu_s) == length(models)
    )
    testthat::expect_true(
      length(test$lambda_m) == length(models)
    )
    testthat::expect_true(
      length(test$mu_m) == length(models)
    )
    testthat::expect_true(
      length(test$lambda_s) == length(models)
    )
    testthat::expect_true(
      length(test$mu_s) == length(models)
    )
    testthat::expect_true(
      length(test$loglik) == length(models)
    )
    testthat::expect_true(
      length(test$df) == length(models)
    )
    testthat::expect_true(
      length(test$conv) == length(models)
    )
    testthat::expect_true(
      length(test$tips_1) == length(models)
    )
    testthat::expect_true(
      length(test$tips_2) == length(models)
    )
    testthat::expect_true(
      length(test$seed) == length(models)
    )
    testthat::expect_true(
      length(unique(test$seed)) == 1
    )
    testthat::expect_true(
      length(test$tips_1) == length(models)
    )
    testthat::expect_true(
      length(test$tips_2) == length(models)
    )
    testthat::expect_true(
      length(test$cond) == length(models)
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
      length(test$optim_lambda_m) == length(models)
    )
    testthat::expect_true(
      length(test$optim_mu_m) == length(models)
    )
    testthat::expect_true(
      length(test$optim_lambda_s) == length(models)
    )
    testthat::expect_true(
      length(test$optim_mu_s) == length(models)
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
      length(test$model) == length(models)
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
    data_file_name <- file.path(
      data_folder,
      paste0(
        "sls_sim",
        "-sim_pars=[0,3-0,2-0,6-0,1]-optim_ids=[1-1-1-1]",
        "-cond=3-n_0=2-ages=[5-2]-seed=",
        seed,
        ".RData"
      )
    )
    testthat::expect_true(
      file.exists(data_file_name)
    )
    suppressWarnings(file.remove(data_file_name))
    # check results file existence
    results_file_name <- file.path(
      results_folder,
      paste0(
        "sls_mle",
        "-sim_pars=[0,3-0,2-0,6-0,1]-optim_ids=[1-1-1-1]",
        "-cond=3-n_0=2-ages=[5-2]-seed=",
        seed,
        ".txt"
      )
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
  # test silent mode and character entry for "models" input
  seed <- 99
  testthat::expect_silent(
    test <- sls_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      l_2 = sim_get_standard_l_2(
        crown_age = 5,
        shift_time = 2
      ),
      start_pars = c(0.2, 0.1, 0.2, 0.1),
      models = "loglik_sls_p",
      verbose = FALSE
    )
  )
  testthat::expect_true(
    is.data.frame(test)
  )
  # check data file existence
  data_file_name <- file.path(
    data_folder,
    paste0(
      "sls_sim",
      "-sim_pars=[0,3-0,2-0,6-0,1]-optim_ids=[1-1-1-1]",
      "-cond=3-n_0=2-ages=[5-2]-seed=",
      seed,
      ".RData"
    )
  )
  testthat::expect_true(
    file.exists(data_file_name)
  )
  suppressWarnings(file.remove(data_file_name))
  # check results file existence
  results_file_name <- file.path(
    results_folder,
    paste0(
      "sls_mle",
      "-sim_pars=[0,3-0,2-0,6-0,1]-optim_ids=[1-1-1-1]",
      "-cond=3-n_0=2-ages=[5-2]-seed=",
      seed,
      ".txt"
    )
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

test_that("it works also for a subset of parameters", {
  seed <- 10
  sim_pars <- c(0.3, 0.2, 0.6, 0.1)
  cond <- 3
  l_2 <- sim_get_standard_l_2(
    crown_age = 5,
    shift_time = 2
  )
  models <- sls_logliks_div()

  test <- sls_main(
    seed = seed,
    sim_pars = sim_pars,
    cond = cond,
    l_2 = l_2,
    start_pars = c(0.2, 0.1, 0.2, 0.1),
    models = models,
    verbose = FALSE,
    optim_ids = c(TRUE, FALSE, FALSE, FALSE)
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
    data_file_name <- file.path(
      data_folder,
      paste0(
        "sls_sim",
        "-sim_pars=[0,3-0,2-0,6-0,1]-optim_ids=[1-0-0-0]",
        "-cond=3-n_0=2-ages=[5-2]-seed=",
        seed,
        ".RData"
      )
    )
    testthat::expect_true(
      file.exists(data_file_name)
    )
    suppressWarnings(file.remove(data_file_name))
    # check results file existence
    results_file_name <- file.path(
      results_folder,
      paste0(
        "sls_mle",
        "-sim_pars=[0,3-0,2-0,6-0,1]-optim_ids=[1-0-0-0]",
        "-cond=3-n_0=2-ages=[5-2]-seed=",
        seed,
        ".txt"
      )
    )
    testthat::expect_true(
      file.exists(results_file_name)
    )
    suppressWarnings(file.remove(results_file_name))
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
      models = "nonsense",
      verbose = FALSE
    ),
    paste0(
      "This is not a likelihood function provided by ",
      get_pkg_name(),
      "!"
    )
  )
})
