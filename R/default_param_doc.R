#' This function does nothing. It is intended to inherit is parameters'
#' documentation.
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param lambdas speciation rates, for all the clades
#' @param mus extinction rate, for all the clades
#' @param Ks carrying capacities, for all the clades
#' @param LS the matrix containing the information about how the subclades are
#' nested into the main clade. See sls_sim.get_standard_LS() for more info.
#' @param LL the collection of all the l tables, for all the clades
#' @param l_matrix_size the initial length of the l matrix. It will be
#' increased if necessary
#' @param data contains all the information about the simulated process
#' @param clade the id of the clade
#' @param delta_n in the Doob-Gillespie algorithm,
#' the change in number of species
#' @param delta_t in the Doob-Gillespie algorithm,
#' the waiting time for the next event to occur
#' @param event the event occurring in the simulated process at a given time
#' @param final_time the final time that you want to consider for the survival
#' of the species considered in the l table
#' @param L the l table
#' @param l_matrix the l table
#' @param pars parameters of the likelihood functions:
#' \itemize{
#'   \item pars[1] is lambda_m, i.e. speciation rate of the main clade;
#'   \item pars[2] is mu_m, i.e. extinction rate of the main clade;
#'   \item pars[3] is lambda_s, i.e. speciation rate of the sub clade;
#'   \item pars[4] is mu_s, i.e. extinction rate of the sub clade;
#' }
#' @param pars_m parameters for the main clade (lambda, mu)
#' @param pars_s parameters for the sub clade (lambda, mu)
#' @param brts branchin times
#' @param brts_m branching times for the Main-clade
#' @param brts_s branching times for the Sub-clade
#' @param n_0 starting number of lineages
#' @param nmax maximum number of lineages to consider
#' @param cond type of conditioning:
#' \itemize{
#'   \item cond = 0 no conditiong;
#'   \item cond = 1 conditions on the survival of crown descendents;
#'   \item cond = 2 not available;
#'   \item cond = 3 conditions on the survival of the subclade and
#'   both crown descendents in the main clade;
#'   \item cond = 4 conditions on the survival of subclade and
#'   on the other crown descendents in the main clade;
#' }
#' @param t time
#' @param ts times
#' @param tbar time left from shift time to the present
#' @param n number of lineages
#' @param fun a function
#' @param k frequencies in the Discrete Fourier Transform (DFT)
#' @param vec a vector or a matrix to be transformed
#' @param E0 starting value for BiSSE's E function
#' @param D0 starting value for BiSSE's D function
#' @param D0s starting values for BiSSE's D functions
#' @param t0 starting time
#' @param tf ending time
#' @param td decoupling time
#' @param tds decoupling times
#' @param LOG set it to TRUE if you desire the output in log form
#' @param lambdaterms set it to TRUE if you desire the powers of lambda
#' in the likelihood
#' @param verbose set it to TRUE if you want to see the outputs on screen
#' @param startpars parameters to start from for the search of the likelihood
#' maximum
#' @param loglik_function the loglik function you want to use
#'
#' @param matrix_size size of the matrix
#' @param lx size of the matrix
#' @param ddep specifies the kind of diversity-dependent model you want to use
#' @param message the message to print

default_params_doc <- function(
  lambda,
  mu,
  lambdas,
  mus,
  Ks,
  LS,
  LL,
  l_matrix_size,
  data,
  clade,
  delta_n,
  delta_t,
  event,
  final_time,
  L,
  l_matrix,
  pars,
  pars_m,
  pars_s,
  brts,
  brts_m,
  brts_s,
  n_0,
  nmax,
  cond,
  t,
  ts,
  tbar,
  n,
  fun,
  k,
  vec,
  E0,
  D0,
  D0s,
  t0,
  tf,
  td,
  tds,
  LOG,
  lambdaterms,
  verbose,
  startpars,
  loglik_function,
  matrix_size,
  lx,
  ddep,
  message
) {
  # Nothing
}
