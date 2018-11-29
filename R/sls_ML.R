#' @title sls Maximum Likelihood
#' @description Calculates ML.
#' @inheritParams default_params_doc
#' @return  best parameters
#' @export
sls_ml <- function(
  loglik_function = sls::loglik_sls_p,
  brts_m,
  brts_s,
  start_pars = c(0.5, 0.3, 0.5, 0.3),
  cond = 3,
  n_0 = 2,
  verbose = TRUE
) {
  if (any(start_pars < 0)) {
    stop("you cannot start from negative parameters")
  }
  failpars <- rep(-1, length(start_pars))
  par_names <- c("lambda_m", "mu_m", "lambda_s", "mu_s")
  out_names <- c(par_names, "loglik", "df", "conv")
  failout  <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
  colnames(failout) <- out_names
  pars <- start_pars

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
    message = paste0("The loglikelihood for the initial parameter values is ", initloglik, "\n"), # nolint
    verbose = verbose
  )
  utils::flush.console()
  if (initloglik == -Inf) {
    cat(
      message = "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n" # nolint
    )
    out2 <- failout
  } else {
    out <- subplex::subplex(
      par = pars2,
      fn = function(x) -fun(pars_transform_back(x))
    ); pars_transform_back(out$par); out[-1]; fun(pars_transform_back(out$par))
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
        df = length(start_pars),
        conv = unlist(out$conv)
      )
      names(out2) <- out_names
    }
  }

  invisible(out2)
}
