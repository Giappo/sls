#' @title Calculates the likelihood for a subset of parameters
#' @description Calculates the likelihood for a subset of parameters
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
sls_loglik_choosepar <- function(
  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  idparsnoshift,
  pars2,
  brts_m,
  brts_s,
  missnumspec = c(0, 0),
  loglik_function = sls::loglik_sls_p
) {
  if (length(missnumspec) != length(c(0, 0))) {
    cat(
      "This likelihood function is meant to work only for missnumspec = c(0,0).\n" # nolint
    )
    loglik <- -Inf
    return(loglik)
  } else {
    if (!all.equal(missnumspec, c(0, 0))) {
      cat(
        "This likelihood function is meant to work only for missnumspec = c(0,0).\n" # nolint
      )
      loglik <- -Inf
      return(loglik)
    }
  }

  namepars <- c("la_M", "mu_M", "K_M", "la_S", "mu_S", "K_S", "t_d")
  Npars <- length(namepars)

  trpars1 <- rep(0, Npars)
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0) {
    trpars1[idparsfix] <- trparsfix
  }
  if (length(idparsnoshift) != 0) {
    trpars1[idparsnoshift] <- trpars1[idparsnoshift - 3]
  }
  brts_m <- -sort(abs(as.numeric(brts_m)), decreasing = TRUE)
  brts_s <- -sort(abs(as.numeric(brts_s)), decreasing = TRUE)

  t_d    <- trpars1[7]
  tsplit <- pars2[4]
  trpars1[7] <- abs(t_d)
  pars2[4]   <- abs(tsplit)

  #Rampal's transformation
  pars1 <- trpars1
  pars1 <- trpars1 / (1 - trpars1)

  names(trpars1) <- namepars
  first_S_branch <- ifelse(length(brts_s) != 0, max(abs(brts_s)), 0)
  if (max(trpars1[1:(Npars - 1)]) > 1 ||
      min(trpars1[1:(Npars - 1)]) < 0 ||
      trpars1[1] <= trpars1[2] ||
      trpars1[4] <= trpars1[5] ||
      abs(pars1[7]) >= abs(pars2[4]) ||
      abs(pars1[7]) <= first_S_branch) {
    loglik <- -Inf
  } else {
    loglik <- loglik_function(
      pars1 = pars1,
      pars2 = pars2,
      brts_m = brts_m,
      brts_s = brts_s,
      missnumspec = missnumspec
    )
    if (is.nan(loglik) || is.na(loglik)) {
      cat("There are parameter values used which cause numerical problems.\n")
      loglik <- -Inf
    }
  }
  return(loglik)
}
