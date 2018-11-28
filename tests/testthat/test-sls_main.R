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
      length(test$tips_m) > 0
    )
    testthat::expect_true(
      length(test$tips_s) > 0
    )
    testthat::expect_true(
      length(test$seed) > 0
    )
    testthat::expect_true(
      length(test$model) > 0
    )
    testthat::expect_true(
      file.exists(
        system.file("extdata", "data", package = "sls")
      )
    )
    testthat::expect_true(
      file.exists(
        system.file(
          "extdata",
          "data",
          paste0("sim_", seed, ".RData"),
          package = "sls"
        )
      )
    )
    testthat::expect_true(
      file.exists(
        system.file(
          "extdata",
          paste0("sls_mle", seed, ".txt"),
          package = "sls"
        )
      )
    )
    testthat::expect_equal(
      read.csv(
        system.file(
          "extdata",
          paste0("sls_mle", seed, ".txt"),
          package = "sls"
        )
      )[, -1],
      test
    )
  }
})
