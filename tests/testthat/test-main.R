context("sim_bd")

test_that("multiplication works", {
  testthat::expect_silent(sim_bd(pars = c(0.6,0.1),time = 10,N0 = 1))
})
