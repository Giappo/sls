context("test_likelihood_formula2")

test_that("multiplication works", {

  if (!is_on_travis()) return()

  Nsims <- 1000000

  sls_data_sets_name <- ls(pattern = "dataset")
  sls_data_sets <- vector("list",length(sls_data_sets_name))
  for (i in 1:length(sls_data_sets_name)){
    sls_data_sets[[i]] <- get(sls_data_sets_name)
    testthat::expect_true(
      test_likelihood_formula2(Nsims = Nsims, dataset = sls_data_sets[[i]])$spread <= 2
    )
  }



})
