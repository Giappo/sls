context("package style")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("package style", {
  if (Sys.info()['sysname'] != "Windows" && is_on_ci()) {
    lintr::expect_lint_free(
      linters = lintr::with_defaults(
        cyclocomp_linter = lintr::cyclocomp_linter(complexity_limit = 30)
      )
    )
  }
})
