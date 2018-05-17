context("test_likelihood_formula2")

test_that("multiplication works", {

  source("R/datasets.R")

  sls_data_sets_name <- ls(pattern = "dataset")
  sls_data_sets <- vector("list",length(sls_data_sets_name))
  for (i in 1:length(sls_data_sets_name)){
    sls_data_sets[[i]] <- get(sls_data_sets_name)
  }

  testthat::expect_true(
    test_likelihood_formula2(Nsims = 10000, dataset = sls_data_sets[[i]])$spread <= 2
  )

})
