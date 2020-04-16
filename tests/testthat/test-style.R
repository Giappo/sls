context("package style")

test_that("package style", {
  lintr::expect_lint_free(
    linters = lintr::with_defaults(
      cyclocomp_linter = lintr::cyclocomp_linter(complexity_limit = 30)
    )
  )
})
