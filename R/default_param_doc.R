#' This function does nothing. It is intended to inherit is parameters'
#' documentation.
#' @param lambdas speciation initiation rates
#' @param mus extinction rates
#' @param ti no idea
#' @param tb no idea
#' @param ts no idea
#' @param tf no idea
#' @param N0 number of lineages at some unknown time
#' @param Nsims number of simulation
#' @param lik_function likelihood function, can be
#'   \code{\link{lik_custom}} or some unknown others
#' @param sim_function simulation function (no idea what that means), can be
#'   \code{\link{sim_custom}} or some unknown others
#' @param input_check if one (why not TRUE?) the input is
#'   checked (why not always?)
#' @note This is an internal function, so it should be marked with
#'   \code{@noRd}. This is not done, as this will disallow all
#'   functions to find the documentation parameters
default_params_doc <- function(
  lambdas,
  mus,
  ti,
  tb,
  ts,
  tf,
  N0,
  Nsims,
  lik_function,
  sim_function,
  input_check
) {
  # Nothing
}
