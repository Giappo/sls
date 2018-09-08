# context("test_likelihood_formula2")
#
# test_that( "Agreement between likelihood and simulations for elementary trees", {
#
#   if (!is_on_travis()) return()
#
#   # Increase for more precision
#   Nsims <- 1
#
#   load_all_data()
#   sls_dataset_shift  <- get("dataset_pure_branching1")
#   sls_dataset_branch <- get("dataset_pure_shifting1")
#
#   testthat::expect_true(
#     test_likelihood_formula(Nsims = Nsims, dataset = sls_dataset_shift)$spread <= 2
#   )
#
#   testthat::expect_true(
#     test_likelihood_formula(Nsims = Nsims, dataset = sls_dataset_branch)$spread <= 2
#   )
#
# })
