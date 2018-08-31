#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
#' @export
sls_ML <- function(loglik_function = sls:::lik_shift_P,
                   initparsopt,
                   parsfix,
                   idparsopt = 1:4,
                   idparsfix = (1:4)[-idparsopt],
                   brtsM,
                   brtsS,
                   cond = 1,
                   optimmethod = 'simplex',
                   pars.transform = 1,
                   res = 10 * (1 + length(c(brtsM, brtsS))),
                   tol = c(1E-3, 1E-4, 1E-6),
                   maxiter = 1000 * round((1.25)^length(idparsopt)),
                   changeloglikifnoconv = FALSE,
                   verbose = TRUE
){# bracket#1
  # - tol = tolerance in optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - maxiter = the maximum number of iterations in the optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - optimmethod = 'simplex' (current default) or 'subplex' (default of previous versions)
  if (!is.numeric(brtsM))
  {
    stop("'brtsM' must be numeric")
  }
  if (!is.numeric(brtsS))
  {
    stop("'brtsS' must be numeric")
  }
  if (length(idparsfix) == 0) {idparsfix <- NULL}
  if (missing(parsfix) && (length(idparsfix) == 0)){parsfix <- idparsfix <- NULL}

  options(warn=-1)
  namepars = c("la1", "mu1", "la2", "mu2"); Npars <- length(namepars); #if you add more parameters to your model just change this
  failpars <- rep(-1, Npars); names(failpars) <- namepars; #those are the parameters that you get if something goes sideways
  if (is.numeric(brtsM) == FALSE)
  {
    cat("The branching times of the main clade should be numeric.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  if (is.numeric(brtsS) == FALSE)
  {
    cat("The branching times of the sub clade should be numeric.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  idpars <- sort(c(idparsopt, idparsfix))
  if ( (sum(idpars == (1:Npars)) != Npars) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)) )
  {
    cat("The parameters to be optimized and/or fixed are incoherent.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(out2)
  }

  if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
  if (verbose == TRUE) {
    cat("You are optimizing",optstr,"\n")
  }
  if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
  if (verbose == TRUE) {
    cat("You are fixing",fixstr,"\n")
    cat("Optimizing the likelihood - this may take a while.","\n")
    flush.console()
  }
  if (pars.transform == 1)
  {
    #Rampal's transformation
    trparsopt = initparsopt/(1 + initparsopt)
    trparsopt[which(initparsopt == Inf)] = 1
    trparsfix = parsfix/(1 + parsfix)
    trparsfix[which(parsfix == Inf)] = 1
  }else
  {
    trparsopt  <- initparsopt
    trparsfix  <- parsfix
  }
  optimpars  <- c(tol, maxiter)
  initloglik <- sls:::sls_loglik_choosepar(loglik_function = loglik_function,
                                           trparsopt = trparsopt, trparsfix = trparsfix,
                                           idparsopt = idparsopt, idparsfix = idparsfix,
                                           brtsM = brtsM, brtsS = brtsS, cond = cond,
                                           pars.transform = pars.transform)
  if (verbose == TRUE) {
    cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
    flush.console()
  }
  if (initloglik == -Inf)
  {# bracket#4
    cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return (invisible(out2))
  }
  if (verbose == TRUE) {
    sink(file = tempfile()) # Sink output here
  }
  out <- DDD:::optimizer(optimmethod = optimmethod, optimpars = optimpars,
                         fun = sls:::sls_loglik_choosepar,
                         loglik_function = loglik_function,
                         trparsopt = trparsopt, trparsfix = trparsfix,
                         idparsopt = idparsopt, idparsfix = idparsfix,
                         brtsM = brtsM, brtsS = brtsS, cond = cond,
                         pars.transform = pars.transform)
  if (verbose == TRUE) {
    sink() # Give back the output
  }
  if (out$conv != 0)
  {# bracket#5
    cat("Optimization has not converged. Try again with different initial values.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  MLtrpars <- as.numeric(unlist(out$par))
  if (pars.transform == 1)
  {
    #Rampal's transformation
    MLpars = MLtrpars/(1 - MLtrpars)
  }else
  {
    MLpars <- MLtrpars
  }
  MLpars1 <- rep(0, Npars); names(MLpars1) <- namepars
  MLpars1[idparsopt] <- MLpars
  if (length(idparsfix) != 0) {MLpars1[idparsfix] <- parsfix}
  ML <- as.numeric(unlist(out$fvalues))
  out2 <- data.frame(t(MLpars1), loglik = ML, df = length(initparsopt), conv = unlist(out$conv))

  tobeprint <- "Maximum likelihood parameter estimates:"
  for (ii in 1:Npars)
  {
    tobeprint <- paste(tobeprint,paste(names(MLpars1[ii]),":",sep = ""),MLpars1[ii])
  }
  if (verbose == TRUE) {
    s1 <- sprintf(tobeprint)
  }
  if(out2$conv != 0 & changeloglikifnoconv == T) {out2$loglik <- -Inf}
  if (verbose == TRUE) {
    s2 <- sprintf('Maximum loglikelihood: %f',ML)
    cat("\n",s1,"\n",s2,"\n\n")
  }

  invisible(out2)
}

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
#' @export
sls_ML2 <- function(loglik_function = sls:::lik_shift_P2,
                    brtsM, brtsS, tsplit,
                    initparsopt = c(0.5, 0.1, 2 * (1 + length(brtsM) + missnumspec[1]), 2 * (1 + length(brtsS) +
                    missnumspec[length(missnumspec)]), (tsplit + max(brtsS))/2),
                    parsfix = NULL,
                    idparsopt = c(1:3, 6:7),
                    idparsfix = NULL,
                    idparsnoshift = (1:7)[c(-idparsopt, (-1)^(length(idparsfix) != 0) * idparsfix)],
                    res = 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec)),
                    ddmodel = 1,
                    missnumspec = c(0,0),
                    cond = 1,
                    soc = 2,
                    tol = c(0.001, 1e-04, 1e-06),
                    maxiter = 1000 * round((1.25)^length(idparsopt)), changeloglikifnoconv = FALSE,
                    optimmethod = "subplex",
                    pars.transform = 1)
{
  namepars <- c("lambda_M", "mu_M", "K_M", "lambda_S", "mu_S", "K_S", "t_d"); Npars <- length(namepars);
  failpars <- rep(-1, Npars); names(failpars) <- namepars; #those are the parameters that you get if something goes sideways
  failout <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  if (missnumspec != c(0,0))
  {
    cat("This likelihood function is meant to work only for missnumspec = c(0,0).\n")
    out2 <- failout
  }

  options(warn = -1)
  brtsM = sort(abs(as.numeric(brtsM)), decreasing = TRUE)
  brtsS = sort(abs(as.numeric(brtsS)), decreasing = TRUE)
  if (cond == 1 & soc == 1)
  {
    cat("Conditioning on survival of a clade with stem age currently not implemented.\n")
    out2 <- failout
  }else
  {
    if (is.numeric(brtsM) == FALSE ||
        is.numeric(brtsS) == FALSE)
    {
      cat("The branching times should be numeric.\n")
      out2 <- failout
    }else
    {
      idparsnoshift = sort(idparsnoshift)
      idpars = sort(c(idparsopt, idparsfix, idparsnoshift))
      if ((prod(idpars == (1:Npars)) != 1) ||
          (length(initparsopt) != length(idparsopt)) ||
          (length(parsfix) != length(idparsfix)))
      {
        cat("The parameters to be optimized, fixed and not shifted are incoherent.\n")
        out2 <- failout
      }else
      {
        # namepars = c("la_M", "mu_M", "K_M", "la_S",
        #              "mu_S", "K_S", "t_d")
        if (length(namepars[idparsopt]) == 0)
        {
          optstr = "nothing"
        }else
        {
          optstr = namepars[idparsopt]
        }
        cat("You are optimizing", optstr, "\n")
        if (length(namepars[idparsfix]) == 0)
        {
          fixstr = "nothing"
        }else
        {
          fixstr = namepars[idparsfix]
        }
        cat("You are fixing", fixstr, "\n")
        if (length(namepars[idparsnoshift]) == 0)
        {
          noshiftstr = "anything"
        }else
        {
          noshiftstr = namepars[idparsnoshift]
        }
        cat("You are not shifting", noshiftstr, "\n")
        cat("Optimizing the likelihood - this may take a while.", "\n")
        flush.console()

        # trparsopt = initparsopt/(1 + initparsopt)
        # trparsopt[which(initparsopt == Inf)] = 1
        # trparsfix = parsfix/(1 + parsfix)
        # trparsfix[which(parsfix == Inf)] = 1
        if (pars.transform == 1)
        {
          #Rampal's transformation
          trparsopt = initparsopt/(1 + initparsopt)
          trparsopt[which(initparsopt == Inf)] = 1
          trparsfix = parsfix/(1 + parsfix)
          trparsfix[which(parsfix == Inf)] = 1
        }else
        {
          trparsopt  <- initparsopt
          trparsfix  <- parsfix
        }

        pars2 = c(res, ddmodel, cond, tsplit, 0, soc, tol, maxiter)
        names(pars2) <- c("res", "ddmodel", "cond", "tsplit", "", "soc", "tol", "maxiter")
        optimpars = c(tol, maxiter)
        initloglik = sls::sls_loglik_choosepar2(loglik_function = loglik_function,
                                                trparsopt = trparsopt,
                                                trparsfix = trparsfix,
                                                idparsopt = idparsopt,
                                                idparsfix = idparsfix,
                                                idparsnoshift = idparsnoshift,
                                                pars2 = pars2, brtsM = brtsM, brtsS = brtsS,
                                                missnumspec = missnumspec,
                                                pars.transform = pars.transform); initloglik
        cat("The loglikelihood for the initial parameter values is", initloglik, "\n")
        flush.console()
        if (initloglik == -Inf)
        {
          cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
          out2 <- failout
        }else
        {
          out = DDD::optimizer(loglik_function = loglik_function,
                               fun = sls::sls_loglik_choosepar2,
                               optimmethod = optimmethod, optimpars = optimpars,
                               trparsopt = trparsopt,
                               trparsfix = trparsfix,
                               idparsopt = idparsopt,
                               idparsfix = idparsfix,
                               idparsnoshift = idparsnoshift,
                               pars2 = pars2, brtsM = brtsM, brtsS = brtsS,
                               missnumspec = missnumspec,
                               pars.transform = pars.transform)
          if (out$conv > 0)
          {
            cat("Optimization has not converged. Try again with different initial values.\n")
            out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = unlist(out$conv))
          }else
          {
            MLtrpars = as.numeric(unlist(out$par))
            if (pars.transform == 1)
            {
              MLpars = MLtrpars/(1 - MLtrpars)
            }else
            {
              MLpars <- MLtrpars
            }
            MLpars1 = rep(0, Npars)
            MLpars1[idparsopt] = MLpars
            if (length(idparsfix) != 0)
            {
              MLpars1[idparsfix] = parsfix
            }
            if (length(idparsnoshift) != 0)
            {
              MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 3]
            }
            if (MLpars1[3] > 10^7)
            {
              MLpars1[3] = Inf
            }
            if (MLpars1[6] > 10^7)
            {
              MLpars1[6] = Inf
            }
            ML = as.numeric(unlist(out$fvalues))
            out2 = data.frame(row.names = "results",
                              lambda_M = MLpars1[1], mu_M = MLpars1[2],
                              K_M = MLpars1[3], lambda_S = MLpars1[4],
                              mu_S = MLpars1[5], K_S = MLpars1[6], t_d = MLpars1[7],
                              loglik = ML, df = length(initparsopt),
                              conv = unlist(out$conv))
            if (out2$conv != 0 & changeloglikifnoconv == TRUE)
            {
              out2$loglik = -Inf
            }
            s1 = sprintf("Maximum likelihood parameter estimates: %f %f %f %f %f %f %f",
                         MLpars1[1], MLpars1[2], MLpars1[3], MLpars1[4],
                         MLpars1[5], MLpars1[6], MLpars1[7])
            s2 = sprintf("Maximum loglikelihood: %f", ML)
            cat("\n", s1, "\n", s2, "\n")
            out$par = list(MLpars1)
          }
        }
      }
    }
  }
  invisible(out2)
}

