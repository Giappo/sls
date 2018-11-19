#' @title Brts to time intervals
#' @author Giovanni Laudanno
#' @description Calculates time intervals given the branching times
#' @inheritParams default_params_doc
#' @return The time intervals
#' @export
brts2time_intervals <- function(
  brts
) {
  time_points <- -unlist(unname(sort(abs(brts), decreasing = TRUE)) )
  branching_times <- -sort(abs(as.numeric(time_points)), decreasing = TRUE)
  births <- c(0, unname(table(branching_times))[-1])
  unique_branching_times <- as.numeric(names(table(branching_times)))
  time_intervals <- c(
    (diff(unique_branching_times)),
    (abs(tail(unique_branching_times, 1)))
  )
  births <- births[-1]
  time_intervals <- c(0, time_intervals)
  births <- c(0, births)
  return(time_intervals)
}

#' @title Transition matrix builder
#' @author Giovanni Laudanno
#' @description Builds the transition matrix to integrate the differential equations of the P-equation
#' @inheritParams default_params_doc
#' @return The transition matrix
#' @export
P_transition_matrix <- function(
  lambda,
  mu,
  matrix_size
) {
  nvec <- 0:(matrix_size - 1)
  M <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  M[row(M) == col(M) + 1] <- M[row(M) == col(M) + 1] +
    lambda * (nvec[1:(matrix_size - 1)])
  M[row(M) == col(M) - 1] <- M[row(M) == col(M) - 1] +
    mu * (nvec[2:matrix_size])
  M[row(M) == col(M)] <- M[row(M) == col(M)] -
    (lambda + mu) * (nvec[1:matrix_size])
  M[matrix_size, matrix_size] <- - M[matrix_size - 1, matrix_size]; M
  testit::assert(colSums(M) < 1e-10)
  return(M)
}
