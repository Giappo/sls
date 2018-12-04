#' @title sls Maximum Likelihood
#' @description Calculates ML.
#' @inheritParams default_params_doc
#' @return  best parameters
#' @export
sls_ml <- function(
  loglik_function = sls::loglik_sls_p,
  brts,
  start_pars = c(0.5, 0.3, 0.5, 0.3),
  n_0 = 2,
  cond = 3,
  verbose = TRUE,
  optim_ids = rep(TRUE, length(start_pars)),
  true_pars = start_pars
) {
  # setup and checks
  par_names <- get_param_names() # nolint internal function
  testit::assert(length(optim_ids) == length(start_pars))
  testit::assert(length(true_pars) == length(start_pars))
  start_pars[!optim_ids] <- true_pars[!optim_ids]
  if (any(start_pars < 0)) {
    stop("You cannot start from negative parameters!")
  }
  out_names <- c(par_names, "loglik", "df", "conv")
  failpars <- rep(-1, length(start_pars))
  failout  <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  colnames(failout) <- out_names

  # define function to optimize
  optim_fun <- function(tr_optim_pars) {
    pars2 <- rep(0, length(start_pars))
    optim_pars <- pars_transform_back(tr_optim_pars) # nolint internal function
    pars2[optim_ids] <- optim_pars
    pars2[!optim_ids] <- true_pars[!optim_ids]

    out <- -loglik_function(
      pars = pars2,
      brts = brts,
      cond = cond,
      n_0 = n_0
    )
    if (verbose == TRUE) {
      printed_values <- paste0(
        c(par_names, "loglik"),
        " = ",
        signif(c(pars2, -out), digits = 5)
      )
      print_this <- paste(printed_values, sep = ",")
      cat(print_this, "\n")
    }
    out
  }

  # initial likelihood
  tr_start_pars <- rep(0, length(start_pars))
  tr_start_pars <- pars_transform_forward(start_pars[optim_ids]) # nolint internal function
  initloglik <- -optim_fun(tr_start_pars)
  utils::flush.console()
  if (initloglik == -Inf) {
    cat(
      message = "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n" # nolint
    )
    out2 <- failout
    return(invisible(out2))
  }

  # maximum likelihood
  out <- subplex::subplex(
    par = tr_start_pars,
    fn = function(x) optim_fun(x)
  )

  # report missed convergence
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
    names(out2) <- out_names
    return(invisible(out2))
  }

  # return mle results
  outpars <- rep(0, length(start_pars))
  outpars[optim_ids] <- pars_transform_back( # nolint internal function
    as.numeric(unlist(out$par))
  )
  outpars[!optim_ids] <- true_pars[!optim_ids]
  names(outpars) <- par_names

  out2 <- data.frame(
    row.names = NULL,
    outpars[1],
    outpars[2],
    outpars[3],
    outpars[4],
    -out$value,
    sum(optim_ids),
    unlist(out$conv)
  )
  names(out2) <- out_names
  return(out2)
}
