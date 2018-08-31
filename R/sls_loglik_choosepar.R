#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
#' @export
sls_loglik_choosepar <- function(loglik_function = sls:::lik_shift_P,
                                 trparsopt,
                                 trparsfix,
                                 idparsopt = 1:4,
                                 idparsfix = (1:4)[-idparsopt],
                                 brtsM,
                                 brtsS,
                                 cond = 1,
                                 pars.transform = 0
                                 # missnumspec = c(0,0)
                                 ){
  #This function provides a likelihood for a subset of parameters. This is built to work inside mbd_minusLL_vs_single_parameter or any optimizer like simplex, optim or subplex
  #idparsopt are the ids of the parameters you want to analyze
  #trparsopt are the values for parameters you want to analyze
  #idparsfix are the ids of the parameters you want to fix
  #trparsfix are the values for parameters you want to fix

  # namepars = c("la_M", "mu_M", "K_M", "la_S", "mu_S", "K_S", "t_d"); Npars <- length(namepars);
  namepars = c("la1", "mu1", "la2", "mu2"); Npars <- length(namepars);
  if (length(trparsopt) == Npars && missing(trparsfix)){trparsfix <- NULL}
  trpars1 = rep(0, Npars)
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0)
  {
    trpars1[idparsfix] <- trparsfix
  }
  if (min(trpars1[1:Npars]) < 0){loglik <- -Inf}else
  {
    if (pars.transform == 1)
    {
      #Rampal's transformation
      pars1 <- trpars1
      # pars1[1:(Npars - 1)] = trpars1[1:(Npars - 1)]/(1 - trpars1[1:(Npars - 1)]) #this for the 7 parameter configuration
      pars1 = trpars1/(1 - trpars1)
    }else
    {
      pars1 <- trpars1
    }
    loglik <- loglik_function(pars = pars1, brtsM = brtsM, brtsS = brtsS, cond = cond)
  }
  if (is.nan(loglik) || is.na(loglik))
  {
    cat("There are parameter values used which cause numerical problems:",trpars1,"\n")
    loglik <- -Inf
  }
  return(loglik)
}


#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
#' @export
sls_loglik_choosepar2 <- function (trparsopt, trparsfix, idparsopt, idparsfix, idparsnoshift,
                                   pars2, brtsM, brtsS, pars.transform = 0, missnumspec = c(0,0),
                                   loglik_function = sls::lik_shift_P2)
{
  methode <- 'analytical'
  if (missnumspec != c(0,0))
  {
    cat("This likelihood function is meant to work only for missnumspec = c(0,0).\n")
    loglik <- -Inf
    return(loglik)
  }
  namepars <- c("la_M", "mu_M", "K_M", "la_S", "mu_S", "K_S", "t_d"); Npars <- length(namepars);

  trpars1 = rep(0, Npars)
  trpars1[idparsopt] = trparsopt
  if (length(idparsfix) != 0)
  {
    trpars1[idparsfix] = trparsfix
  }
  if (length(idparsnoshift) != 0)
  {
    trpars1[idparsnoshift] = trpars1[idparsnoshift - 3]
  }
  brtsM = -sort(abs(as.numeric(brtsM)), decreasing = TRUE)
  brtsS = -sort(abs(as.numeric(brtsS)), decreasing = TRUE)

  t_d    <- trpars1[7]
  tsplit <- pars2[4]
  t_d    <- abs(t_d)    * sign(brtsM[1])
  tsplit <- abs(tsplit) * sign(brtsM[1])
  # pars1[7]   <- t_d
  trpars1[7] <- t_d
  pars2[4]   <- tsplit

  if (pars.transform == 1)
  {
    #Rampal's transformation
    pars1 <- trpars1
    # pars1[1:(Npars - 1)] = trpars1[1:(Npars - 1)]/(1 - trpars1[1:(Npars - 1)]) #this for the 7 parameter configuration. If you change this check also sls_ML2
    pars1 = trpars1/(1 - trpars1)
  }else
  {
    pars1 <- trpars1
  }

  if (max(trpars1[1:(Npars - 1)]) > 1 ||
      min(trpars1[1:(Npars - 1)]) < 0 ||
      trpars1[1] <= trpars1[2] ||
      trpars1[4] <= trpars1[5] ||
      abs(pars1[7]) >= abs(pars2[4]) ||
      abs(pars1[7]) <= min(brtsS))
  {
    loglik = -Inf
  }else
  {
    loglik = loglik_function(pars1 = pars1, pars2 = pars2,
                brtsM = brtsM, brtsS = brtsS,
                missnumspec = missnumspec)
    if (is.nan(loglik) || is.na(loglik))
    {
      cat("There are parameter values used which cause numerical problems.\n")
      loglik = -Inf
    }
  }
  return(loglik)
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
#' @export
dd_KI_loglik_choosepar2 <- function (trparsopt, trparsfix, idparsopt, idparsfix, idparsnoshift,
                                     pars2, brtsM, brtsS, pars.transform = 0, missnumspec = c(0,0),
                                     fun = sls::lik_shift_DDD2)
{
  methode <- 'analytical'
  if (missnumspec != c(0,0))
  {
    cat("This likelihood function is meant to work only for missnumspec = c(0,0).\n")
    loglik <- -Inf
    return(loglik)
  }
  namepars <- c("la_M", "mu_M", "K_M", "la_S", "mu_S", "K_S", "t_d"); Npars <- length(namepars);

  trpars1 = rep(0, Npars)
  trpars1[idparsopt] = trparsopt
  if (length(idparsfix) != 0)
  {
    trpars1[idparsfix] = trparsfix
  }
  if (length(idparsnoshift) != 0)
  {
    trpars1[idparsnoshift] = trpars1[idparsnoshift - 3]
  }
  brtsM = -sort(abs(as.numeric(brtsM)), decreasing = TRUE)
  brtsS = -sort(abs(as.numeric(brtsS)), decreasing = TRUE)

  # t_d    <- trpars1[7]
  # tsplit <- pars2[4]
  # t_d    <- abs(t_d)    * sign(brtsM[1])
  # tsplit <- abs(tsplit) * sign(brtsM[1])
  # # pars1[7]   <- t_d
  # trpars1[7] <- t_d
  # pars2[4]   <- tsplit

  if (pars.transform == 1)
  {
    #Rampal's transformation
    pars1 <- trpars1
    # pars1[1:(Npars - 1)] = trpars1[1:(Npars - 1)]/(1 - trpars1[1:(Npars - 1)]) #this for the 7 parameter configuration. If you change this check also sls_ML2
    pars1 = trpars1/(1 - trpars1)
  }else
  {
    pars1 <- trpars1
  }

  if (max(trpars1[1:(Npars - 1)]) > 1 ||
      min(trpars1[1:(Npars - 1)]) < 0 ||
      trpars1[1] <= trpars1[2] ||
      trpars1[4] <= trpars1[5] ||
      abs(pars1[7]) >= abs(pars2[4]) ||
      abs(pars1[7]) <= min(brtsS))
  {
    loglik = -Inf
  }else
  {
    loglik = fun(pars1 = pars1, pars2 = pars2,
                 brtsM = brtsM, brtsS = brtsS,
                 missnumspec = missnumspec)
    if (is.nan(loglik) || is.na(loglik))
    {
      cat("There are parameter values used which cause numerical problems.\n")
      loglik = -Inf
    }
  }
  return(loglik)
}
