# context("arrange_times_matrix & times_matrix2t_coordinates")
#
# test_that( "Test 'time matrices' and 'time points' conversion functions", {
#
#   # if (!is_on_travis()) return()
#
#   load_all_data(the.environment = environment())
#   data.sets <- ls(pattern = "dataset_",envir = environment())
#
#   for (i in 1:length(data.sets))
#   {
#     data         <- get(data.sets[[i]])
#     # original_coords$lambdas <- NULL
#     # original_coords$mus     <- NULL
#
#     # time_matrix     <- arrange_times_matrix(ti = original_coords$ti, tb = original_coords$tb,
#                                             # ts = original_coords$ts, tf = original_coords$tf)
#     times_matrix    <- data$times_matrix
#     tcoords         <- times_matrix2t_coordinates(times_matrix = times_matrix)
#
#     testthat::expect_equal(
#       unname(tcoords$ti), unname(original_coords$ti)
#     )
#     testthat::expect_equal(
#       unname(tcoords$tb), unname(original_coords$tb)
#     )
#     testthat::expect_equal(
#       unname(tcoords$ts), unname(original_coords$ts)
#     )
#     testthat::expect_equal(
#       unname(tcoords$tf), unname(original_coords$tf)
#     )
#   }
#
# })

# context("lik_custom & lik_custom_split")
#
# test_that( "Agreement between lik_custom and lik_custom_split for elementary trees", {
#
#   #they are not equal! Fortunately!
#
#   # if (!is_on_travis()) return()
#
#   # load_all_data(the.environment = environment())
#   # data.sets <- ls(pattern = "dataset_",envir = environment())
#   #
#   # for (i in 1:length(data.sets))
#   # {
#   #   data <- get(data.sets[[i]])
#   #
#   #   out1 <- lik_custom(lambdas = data$lambdas, mus = data$mus, ti = data$ti, tb = data$tb, ts = data$ts, tf = data$tf)
#   #   out2 <- lik_custom_split(lambdas = data$lambdas, mus = data$mus, ti = data$ti, tb = data$tb, ts = data$ts, tf = data$tf)
#   #   testthat::expect_equal(
#   #     out1, out2
#   #   )
#   # }
#
# })

# context("lik_custom_split & lik_custom_split2")
#
# test_that( "Agreement between lik_custom_split and lik_custom_split2", {
#
#   # if (!is_on_travis()) return()
#
#   load_all_data(the.environment = environment())
#   data.sets <- ls(pattern = "dataset_",envir = environment())
#
#   for (i in 1:length(data.sets))
#   {
#     dataset <- get(data.sets[[i]])
#
#     out1 <- lik_custom_split(dataset = dataset)
#     out2 <- lik_custom_split2(dataset = dataset)
#
#     testthat::expect_equal(
#       out1, out2
#     )
#
#     # out3 <- lik_custom_split3(dataset = dataset)
#     # testthat::expect_equal(
#     # out2, out3
#     # )
#   }
#
# })

# context("test_likelihood_formula2")
#
# test_that( "Agreement between likelihood and simulations for elementary trees", {
#
#   if (!is_on_travis()) return()
#
#   Nsims <- 1000000
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
