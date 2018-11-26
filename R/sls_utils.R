#' @title Transition matrix builder
#' @author Giovanni Laudanno
#' @description Builds the transition matrix to integrate the differential equations of the P-equation
#' @inheritParams default_params_doc
#' @return The transition matrix
#' @export
p_transition_matrix <- function(
  lambda,
  mu,
  matrix_size
) {
  nvec <- 0:(matrix_size - 1)
  m <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  m[row(m) == col(m) + 1] <- m[row(m) == col(m) + 1] +
    lambda * (nvec[1:(matrix_size - 1)])
  m[row(m) == col(m) - 1] <- m[row(m) == col(m) - 1] +
    mu * (nvec[2:matrix_size])
  m[row(m) == col(m)] <- m[row(m) == col(m)] -
    (lambda + mu) * (nvec[1:matrix_size])
  m[matrix_size, matrix_size] <- - m[matrix_size - 1, matrix_size]; m
  testit::assert(colSums(m) < 1e-10)
  return(m)
}

#' @title cat2
#' @author Giovanni Laudanno
#' @description If verbose == TRUE cats the message, otherwise stays silent
#' @inheritParams default_params_doc
#' @return prints on screen
#' @export
cat2 <- function(
  message,
  verbose
) {
 if (verbose == TRUE) {
   cat(message)
 } else {
   return()
 }
}

#' @title Check input
#' @author Giovanni Laudanno
#' @description It checks the inputs
#' @inheritParams default_params_doc
#' @return nothing
#' @export
sls_check_input <- function(
  brts_m,
  brts_s,
  cond,
  n_0,
  n_max
) {
  if (length(brts_m) <= 0) {
    stop("main clade branching times cannot be an empty vector")
  }
  if (length(brts_s) <= 0) {
    stop("sub clade branching times cannot be an empty vector")
  }
  if (any(brts_m < 0)) {
    stop("all the branching times for the main clade have to be non negative")
  }
  if (any(brts_s < 0)) {
    stop("all the branching times for the sub clade have to be non negative")
  }
  if (!(cond %in% sls_conds())) {
   stop("this conditioning is not implemented")
  }
  if (!(n_0 %in% sls_n_0s())) {
   stop("this n_0 is not implemented")
  }
  if (n_max <= 0) {
   stop("it's not going to work with maximum species set to 0 or less")
  }
}

#' @title Conditionings
#' @author Giovanni Laudanno
#' @description Gives the conditionings accepted by sls
#' @inheritParams default_params_doc
#' @return the conditionings
#' @export
sls_conds <- function() {
  conds <- c(3, 4)
  conds
}

#' @title Starting species
#' @author Giovanni Laudanno
#' @description Gives the amount of starting species accepted by sls
#' @inheritParams default_params_doc
#' @return the possible n_0s
#' @export
sls_n_0s <- function() {
  n_0s <- c(2)
  n_0s
}

#' @title Logliks with division
#' @author Giovanni Laudanno
#' @description Get all the loglik functions with division
#' @inheritParams default_params_doc
#' @return loglik functions with division in sls
#' @export
sls_logliks_div <- function() {
  fun_list <- ls("package:sls")
  div_funs <- fun_list[sapply(
    fun_list, function(x)
      any(grepl("loglik_sls", x)) &
      !any(grepl("nodiv", x))
  )]
  div_funs
}

#' @title Logliks with no division
#' @author Giovanni Laudanno
#' @description Get all the loglik functions with no division
#' @inheritParams default_params_doc
#' @return loglik functions with no division in sls
#' @export
sls_logliks_nodiv <- function() {
  fun_list <- ls("package:sls")
  nodiv_funs <- fun_list[sapply(
    fun_list, function(x)
      any(grepl("loglik", x)) &
      any(grepl("nodiv", x))
  )]
  nodiv_funs
}
