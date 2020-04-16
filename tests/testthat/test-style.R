context("package style")

test_that("package style", {
  lintr::expect_lint_free(
    linters = lintr::cyclocomp_linter(complexity_limit = 30)
  )
})
