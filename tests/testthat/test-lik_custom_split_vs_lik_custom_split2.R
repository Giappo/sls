context("lik_custom_split & lik_custom_split2")

test_that( "Agreement between lik_custom_split and lik_custom_split2", {

  # if (!is_on_travis()) return()

  load_all_data(the.environment = environment())
  data.sets <- ls(pattern = "dataset_",envir = environment())

  for (i in 1:length(data.sets))
  {
    dataset <- get(data.sets[[i]])

    out1 <- lik_custom_split(dataset = dataset)
    out2 <- lik_custom_split2(dataset = dataset)
    # out3 <- lik_custom_split3(lambdas = data$lambdas, mus = data$mus, ti = data$ti, tb = data$tb, ts = data$ts, tf = data$tf)
    testthat::expect_equal(
      out1, out2
    )
    # testthat::expect_equal(
      # out2, out3
    # )
  }

})
