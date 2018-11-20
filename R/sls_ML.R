#' @title sls Maximum Likelihood
#' @description Calculates ML.
#' @inheritParams default_params_doc
#' @return  best parameters
#' @export
sls_ml <- function(
  loglik_function = sls::loglik_slsP,
  brts_m,
  brts_s,
  startpars = c(0.5, 0.3, 0.5, 0.3),
  cond = 3,
  n_0 = 2,
  verbose = 1
) {
  if (any(startpars < 0)) {
    stop("you cannot start from negative parameters")
  }
  failpars <- rep(-1, length(startpars))
  par_names <- c("lambda_s", "mu_m", "lambda_s", "mu_s")
  failout  <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  pars <- startpars

  #Rampal's transformation
  pars2 <- pars / (1 + pars)
  pars2[which(pars == Inf)] <- 1

  fun <- function(pars) {
    loglik_function(
      pars_m = pars[1:2],
      pars_s = pars[3:4],
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_0 = n_0
    )
  }

  initloglik <- fun(pars); pars; initloglik
  cat2(
    message = paste0("The loglikelihood for the initial parameter values is", initloglik, "\n"), # nolint
    verbose = verbose
  )
  utils::flush.console()
  if (initloglik == -Inf) {
    cat2(
      message = paste0("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n"), # nolint
      verbose = verbose
    )
    out2 <- failout
  } else {
    out <- subplex::subplex(
      par = pars2,
      fn = function(x) -fun(x)
    ); out; fun(out$par)
    if (out$conv > 0) {
      cat2(
        "Optimization has not converged. Try again with different initial values.\n", # nolint
        verbose = verbose
      )
      out2 <- data.frame(
        t(failpars),
        loglik = -1,
        df = -1,
        conv = unlist(out$conv)
      )
    } else {
      outpars2 <- as.numeric(unlist(out$par))
      outpars <- outpars2 / (1 - outpars2)
      names(outpars) <- par_names

      out2 <- data.frame(
        row.names = "results",
        lambda_m = outpars[1],
        mu_m = outpars[2],
        lambda_s = outpars[3],
        mu_s = outpars[4],
        loglik = out$value,
        df = length(startpars),
        conv = unlist(out$conv)
      )
    }
  }

  invisible(out2)
}

  #' #' @title sls Maximum Likelihood
  #' #' @description Calculates ML.
  #' #' @details best parameters
  #' #' @export
  #' sls_ML <- function(loglik_function = sls::loglik_slsP,
  #'                     brts_m, brts_s, tsplit,
  #'                     initparsopt = c(0.5, 0.1, 2 * (1 + length(brts_m) + missnumspec[1]), 2 * (1 + length(brts_s) +
  #'                     missnumspec[length(missnumspec)]), (tsplit + max(brts_s))/2),
  #'                     parsfix = NULL,
  #'                     idparsopt = c(1:3, 6:7),
  #'                     idparsfix = NULL,
  #'                     idparsnoshift = (1:7)[c(-idparsopt, (-1)^(length(idparsfix) != 0) * idparsfix)],
  #'                     res = 10 * (1 + length(c(brts_m, brts_s)) + sum(missnumspec)),
  #'                     ddmodel = 1,
  #'                     missnumspec = c(0,0),
  #'                     cond = 1,
  #'                     soc = 2,
  #'                     tol = c(0.001, 1e-04, 1e-06),
  #'                     maxiter = 1000 * round((1.25)^length(idparsopt)),
  #'                     changeloglikifnoconv = FALSE,
  #'                     optimmethod = 'subplex')
  #' {
  #'   namepars <- c("lambda_M", "mu_M", "K_M", "lambda_S", "mu_S", "K_S", "t_d"); Npars <- length(namepars);
  #'   failpars <- rep(-1, Npars); names(failpars) <- namepars; #those are the parameters that you get if something goes sideways
  #'   failout  <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  #'   if (length(missnumspec) != length(c(0,0)))
  #'   {
  #'     cat("This likelihood function is meant to work only for missnumspec = c(0,0).\n")
  #'     out2 <- failout
  #'   }else
  #'   {
  #'     if (!all.equal(missnumspec, c(0,0)))
  #'     {
  #'       cat("This likelihood function is meant to work only for missnumspec = c(0,0).\n")
  #'       out2 <- failout
  #'     }
  #'   }
  #'
  #'   options(warn = -1)
  #'   brts_m = sort(abs(as.numeric(brts_m)), decreasing = TRUE)
  #'   brts_s = sort(abs(as.numeric(brts_s)), decreasing = TRUE)
  #'   if (cond != 0 & soc == 1)
  #'   {
  #'     cat("Conditioning on survival of a clade with stem age currently not implemented.\n")
  #'     out2 <- failout
  #'   }else
  #'   {
  #'     if (is.numeric(brts_m) == FALSE ||
  #'         is.numeric(brts_s) == FALSE)
  #'     {
  #'       cat("The branching times should be numeric.\n")
  #'       out2 <- failout
  #'     }else
  #'     {
  #'       idparsnoshift = sort(idparsnoshift)
  #'       idpars = sort(c(idparsopt, idparsfix, idparsnoshift))
  #'       if ((prod(idpars == (1:Npars)) != 1) ||
  #'           (length(initparsopt) != length(idparsopt)) ||
  #'           (length(parsfix) != length(idparsfix)))
  #'       {
  #'         cat("The parameters to be optimized, fixed and not shifted are incoherent.\n")
  #'         out2 <- failout
  #'       }else
  #'       {
  #'         # namepars = c("la_M", "mu_M", "K_M", "la_S",
  #'         #              "mu_S", "K_S", "t_d")
  #'         if (length(namepars[idparsopt]) == 0)
  #'         {
  #'           optstr = "nothing"
  #'         }else
  #'         {
  #'           optstr = namepars[idparsopt]
  #'         }
  #'         cat("You are optimizing", optstr, "\n")
  #'         if (length(namepars[idparsfix]) == 0)
  #'         {
  #'           fixstr = "nothing"
  #'         }else
  #'         {
  #'           fixstr = namepars[idparsfix]
  #'         }
  #'         cat("You are fixing", fixstr, "\n")
  #'         if (length(namepars[idparsnoshift]) == 0)
  #'         {
  #'           noshiftstr = "anything"
  #'         }else
  #'         {
  #'           noshiftstr = namepars[idparsnoshift]
  #'         }
  #'         cat("You are not shifting", noshiftstr, "\n")
  #'         cat("Optimizing the likelihood - this may take a while.", "\n")
  #'         utils::flush.console()
  #'
  #'         #Rampal's transformation
  #'         trparsopt = initparsopt/(1 + initparsopt)
  #'         trparsopt[which(initparsopt == Inf)] = 1
  #'         trparsfix = parsfix/(1 + parsfix)
  #'         trparsfix[which(parsfix == Inf)] = 1
  #'
  #'         pars2 = c(res, ddmodel, cond, tsplit, 0, soc, tol, maxiter)
  #'         names(pars2) <- c("res", "ddmodel", "cond", "tsplit", "", "soc", "tol", "maxiter")
  #'         optimpars = c(tol, maxiter)
  #'         initloglik = sls::sls_loglik_choosepar(loglik_function = loglik_function,
  #'                                                trparsopt = trparsopt,
  #'                                                trparsfix = trparsfix,
  #'                                                idparsopt = idparsopt,
  #'                                                idparsfix = idparsfix,
  #'                                                idparsnoshift = idparsnoshift,
  #'                                                pars2 = pars2, brts_m = brts_m, brts_s = brts_s,
  #'                                                missnumspec = missnumspec); initloglik
  #'         cat("The loglikelihood for the initial parameter values is", initloglik, "\n")
  #'         utils::flush.console()
  #'         if (initloglik == -Inf)
  #'         {
  #'           cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
  #'           out2 <- failout
  #'         }else
  #'         {
  #'           out = DDD::optimizer(loglik_function = loglik_function,
  #'                                fun = sls::sls_loglik_choosepar,
  #'                                optimmethod = optimmethod,
  #'                                optimpars = optimpars,
  #'                                trparsopt = trparsopt,
  #'                                trparsfix = trparsfix,
  #'                                idparsopt = idparsopt,
  #'                                idparsfix = idparsfix,
  #'                                idparsnoshift = idparsnoshift,
  #'                                pars2 = pars2, brts_m = brts_m, brts_s = brts_s,
  #'                                missnumspec = missnumspec)
  #'           if (out$conv > 0)
  #'           {
  #'             cat("Optimization has not converged. Try again with different initial values.\n")
  #'             out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = unlist(out$conv))
  #'           }else
  #'           {
  #'             MLtrpars = as.numeric(unlist(out$par))
  #'             MLpars = MLtrpars/(1 - MLtrpars)
  #'             MLpars1 = rep(0, Npars)
  #'             MLpars1[idparsopt] = MLpars
  #'             if (length(idparsfix) != 0)
  #'             {
  #'               MLpars1[idparsfix] = parsfix
  #'             }
  #'             if (length(idparsnoshift) != 0)
  #'             {
  #'               MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 3]
  #'             }
  #'             if (MLpars1[3] > 10^7)
  #'             {
  #'               MLpars1[3] = Inf
  #'             }
  #'             if (MLpars1[6] > 10^7)
  #'             {
  #'               MLpars1[6] = Inf
  #'             }
  #'             ML = as.numeric(unlist(out$fvalues))
  #'             out2 = data.frame(row.names = "results",
  #'                               lambda_M = MLpars1[1], mu_M = MLpars1[2],
  #'                               K_M = MLpars1[3], lambda_S = MLpars1[4],
  #'                               mu_S = MLpars1[5], K_S = MLpars1[6], t_d = MLpars1[7],
  #'                               loglik = ML, df = length(initparsopt),
  #'                               conv = unlist(out$conv))
  #'             if (out2$conv != 0 & changeloglikifnoconv == TRUE)
  #'             {
  #'               out2$loglik = -Inf
  #'             }
  #'             s1 = sprintf("Maximum likelihood parameter estimates: %f %f %f %f %f %f %f",
  #'                          MLpars1[1], MLpars1[2], MLpars1[3], MLpars1[4],
  #'                          MLpars1[5], MLpars1[6], MLpars1[7])
  #'             s2 = sprintf("Maximum loglikelihood: %f", ML)
  #'             cat("\n", s1, "\n", s2, "\n")
  #'             out$par = list(MLpars1)
  #'           }
  #'         }
  #'       }
  #'     }
  #'   }
  #'   invisible(out2)
  #' }
