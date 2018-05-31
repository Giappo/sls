context("lik_custom_split2 & lik_custom_split3")

test_that( "Agreement between lik_custom_split2 and lik_custom_split3", {

  # if (!is_on_travis()) return()

  # load_all_data(the.environment = environment())
  # data.sets <- ls(pattern = "dataset_",envir = environment())
  #
  # for (i in 1:length(data.sets))
  # {
  #   data <- get(data.sets[[i]])
  #
  #   # out1 <- lik_custom_split(lambdas = data$lambdas, mus = data$mus, ti = data$ti, tb = data$tb, ts = data$ts, tf = data$tf)
  #   out2 <- lik_custom_split2(lambdas = data$lambdas, mus = data$mus, ti = data$ti, tb = data$tb, ts = data$ts, tf = data$tf)
  #   out3 <- lik_custom_split3(lambdas = data$lambdas, mus = data$mus, ti = data$ti, tb = data$tb, ts = data$ts, tf = data$tf)
  #   # testthat::expect_equal(
  #     # out1, out2
  #   # )
  #   testthat::expect_equal(
  #     out2, out3
  #   )
  # }

})
