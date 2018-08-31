#' This function does nothing. It is intended to inherit is parameters'
#' documentation.
#' @param pars parameters of the likelihood functions:
#' \itemize{
#'   \item pars[1] is lambdaM, i.e. speciation rate of the main clade;
#'   \item pars[2] is muM, i.e. extinction rate of the main clade;
#'   \item pars[3] is lambdaS, i.e. speciation rate of the sub clade;
#'   \item pars[4] is muS, i.e. extinction rate of the sub clade;
#' }
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param brts branchin times
#' @param brtsM branching times for the Main-clade
#' @param brtsS branching times for the Sub-clade
#' @param t time
#' @param tbar time left from shift time to the present
#' @param n number of lineages
#' @param nmax maximum number of lineages to consider
#' @param matrix_size size of the matrix
#' @param lx size of the matrix
#' @param ddep specifies the kind of diversity-dependent model you want to use
#' @param k frequencies in discrete fourier transform
#' @param cond type of conditioning:
#' \itemize{
#'   \item cond = 0 no conditiong;
#'   \item cond = 1 conditions on the survival of crown descendents;
#'   \item cond = 2 conditions on the survival of subclade and on the other crown descendents in the main clade;
#'   \item cond = 3 conditions on the survival of the subclade and both crown descendents in the main clade;
#' }
#' @param vec a vector or a matrix to be transformed
default_params_doc <- function(
  lambda,
  mu,
  t,
  n,
  tbar,
  pars,
  brts,
  brtsM,
  brtsS,
  nmax,
  matrix_size,
  lx,
  ddep,
  k,
  cond,
  vec
) {
  # Nothing
}
