context("arrange_times_matrix & times_matrix2t_coordinates")

test_that( "Test 'time matrices' and 'time points' conversion functions", {

  # if (!is_on_travis()) return()

  load_all_data(the.environment = environment())
  data.sets <- ls(pattern = "dataset_",envir = environment())

  for (i in 1:length(data.sets))
  {
    original_coords         <- get(data.sets[[i]])
    original_coords$lambdas <- NULL
    original_coords$mus     <- NULL

    time_matrix     <- arrange_times_matrix(ti = original_coords$ti, tb = original_coords$tb,
                                            ts = original_coords$ts, tf = original_coords$tf)
    tcoords         <- times_matrix2t_coordinates(times_matrix = time_matrix)

    testthat::expect_equal(
      unname(tcoords$ti), unname(original_coords$ti)
    )
    testthat::expect_equal(
      unname(tcoords$tb), unname(original_coords$tb)
    )
    testthat::expect_equal(
      unname(tcoords$ts), unname(original_coords$ts)
    )
    testthat::expect_equal(
      unname(tcoords$tf), unname(original_coords$tf)
    )
  }

})
