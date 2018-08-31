context("lik_custom_split & lik_custom_split2")

test_that( "Agreement between lik_custom_split and lik_custom_split2", {

  load_all_data(the.environment = environment())
  data.sets <- ls(pattern = "dataset_",envir = environment())

  for (i in 1:length(data.sets))
  {
    dataset <- get(data.sets[[i]])

    out1 <- lik_custom_split(dataset = dataset)
    out2 <- lik_custom_split2(dataset = dataset)

    testthat::expect_equal(
      out1, out2
    )

    # out3 <- lik_custom_split3(dataset = dataset)
    # testthat::expect_equal(
    # out2, out3
    # )
  }

})
