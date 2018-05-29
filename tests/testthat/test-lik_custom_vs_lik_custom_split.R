context("test_likelihood_formula2")

test_that( "Agreement between lik_custom and lik_custom_split for elementary trees", {

  # if (!is_on_travis()) return()

  load_all_data()
  data.sets <- ls(pattern = "dataset_")

  for (i in 1:length(data.sets))
  {
  data <- get(data.sets[[i]])

  out1 <- lik_custom(lambdas = data$lambdas, mus = data$mus, ti = data$ti, tb = data$tb, ts = data$ts, tf = data$tf)
  out2 <- lik_custom_split(lambdas = data$lambdas, mus = data$mus, ti = data$ti, tb = data$tb, ts = data$ts, tf = data$tf)
  testthat::expect_equal(
    out1, out2
  )
  }

})
