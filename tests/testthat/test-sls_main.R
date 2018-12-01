context("sls_main")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("use", {
  seed_interval <- 6:(6 + 9 * is_on_ci()) # 6 is critical
  for (seed in seed_interval) {
    sim_pars <- c(0.3, 0.2, 0.6, 0.1)
    cond <- 3
    test <- sls_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      l_2 = sim_get_standard_l_2(
        crown_age = 5,
        shift_time = 2
      ),
      start_pars = c(0.2, 0.1, 0.2, 0.1),
      models = sls_logliks_div(),
      verbose = FALSE
    )
    testthat::expect_true(
      is.data.frame(test)
    )
    testthat::expect_true(
      length(test$sim_lambda_m) > 0
    )
    testthat::expect_true(
      length(test$sim_mu_m) > 0
    )
    testthat::expect_true(
      length(test$sim_lambda_s) > 0
    )
    testthat::expect_true(
      length(test$sim_mu_s) > 0
    )
    testthat::expect_true(
      length(test$lambda_m) > 0
    )
    testthat::expect_true(
      length(test$mu_m) > 0
    )
    testthat::expect_true(
      length(test$lambda_s) > 0
    )
    testthat::expect_true(
      length(test$mu_s) > 0
    )
    testthat::expect_true(
      length(test$loglik) > 0
    )
    testthat::expect_true(
      length(test$df) > 0
    )
    testthat::expect_true(
      length(test$conv) > 0
    )
    testthat::expect_true(
      length(test$tips_1) > 0
    )
    testthat::expect_true(
      length(test$tips_2) > 0
    )
    testthat::expect_true(
      length(test$seed) > 0
    )
    testthat::expect_true(
      length(test$model) > 0
    )

    pkg_name <- get_pkg_name() # nolint internal function
    # test file saving
    if (.Platform$OS.type == "windows") {
      sim_path  <- system.file("extdata", package = pkg_name)
    } else {
      sim_path  <- getwd()
    }
    # check data_path folder existence
    data_path <- file.path(sim_path, "data")
    testthat::expect_true(
      file.exists(data_path)
    )
    # check data file existence
    data_file_name <- file.path(
      data_path,
      paste0(pkg_name, "_sim_", seed, ".RData")
    )
    testthat::expect_true(
      file.exists(data_file_name)
    )
    # check results file existence
    results_file_name <- file.path(
      sim_path,
      paste0(pkg_name, "_mle_", seed, ".txt")
    )
    testthat::expect_true(
      file.exists(results_file_name)
    )
    # check if saved results are the right ones
    testthat::expect_equal(
      utils::read.csv(results_file_name)[, -1],
      test
    )
  }
  # test silent mode and character entry for "models" input
  testthat::expect_silent(
    out <- sls_main(
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
    is.data.frame(out)
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
