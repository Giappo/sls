#' @title P-likelihood dummy
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convolving Nee's functions
#'  with lambda_s == lambda_m and mu_s == mu_m
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_sls_p_dummy <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  n_max = 1e2
) {
  pars2 <- c(pars, pars)
  loglik <- loglik_sls_p(
    pars = pars2,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    n_max = n_max
  )
  loglik
}

#' @title P-likelihood (with no division)
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convoluting Nee's functions
#'  with lambda_s == lambda_m and mu_s == mu_m.
#'  There is no division. It should yield the same likelihood as DDD.
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_sls_p_dummy_nodiv <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  n_max = 1e2
) {
  pars2 <- c(pars, pars)
  loglik <- loglik_sls_p_nodiv(
    pars = pars2,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    n_max = n_max
  )
  loglik
}

#' @title Logliks with dummy shift
#' @author Giovanni Laudanno
#' @description Get all the loglik functions with dummy shift
#' @inheritParams default_params_doc
#' @return loglik functions with dummy shift
#' @export
sls_logliks_dummy <- function() {
  dummy_funs <- c(loglik_sls_p_dummy, loglik_sls_p_dummy_nodiv)
  dummy_funs
}

#' @title sls Maximum Likelihood
#' @description Calculates ML.
#' @inheritParams default_params_doc
#' @return  best parameters
#' @export
sls_ml_dummy <- function(
  loglik_function = sls::loglik_sls_p_dummy,
  brts,
  start_pars = c(0.5, 0.3),
  n_0 = 2,
  cond = 3,
  optim_ids = rep(TRUE, length(start_pars)),
  true_pars = start_pars,
  verbose = TRUE
) {
  # setup and checks
  par_names <- get_param_names()[seq_along(start_pars)] # nolint internal function
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
  outpars_s <- outpars
  names(outpars_s) <- c("lambda_s", "mu_s")

  out2 <- data.frame(
    row.names = NULL,
    outpars[1],
    outpars[2],
    outpars_s[1],
    outpars_s[2],
    -out$value,
    sum(optim_ids),
    unlist(out$conv)
  )
  names(out2) <- c(
    out_names[1:2],
    names(outpars_s),
    out_names[3:length(out_names)]
  )
  return(out2)
}

#' @title Run the full sls routine with dummy shift
#' @description Simulate an sls process and infer parameters maximizing the
#' likelihood(s) for given function(s) with dummy shift
#' @inheritParams default_params_doc
#' @details mle inference with dummy shift
#' @export
sls_main_dummy <- function(
  sim_pars,
  cond = 3,
  l_2 = sim_get_standard_l_2(
    crown_age = 5,
    shift_time = 2
  ),
  seed,
  start_pars = c(0.2, 0.1),
  optim_ids = rep(TRUE, length(start_pars)),
  loglik_functions = sls_logliks_dummy(),
  project_folder = NULL,
  verbose = FALSE
) {
  # check formats
  sim_pars <- as.numeric(sim_pars)
  start_pars <- as.numeric(start_pars)
  cond <- as.numeric(cond)
  seed <- as.numeric(seed)
  names(sim_pars) <- sls::get_param_names()

  # specific set up
  lambdas <- sim_pars[c(1, 3)]
  mus <- sim_pars[c(2, 4)]
  ks <- c(Inf, Inf)
  t_0s <- l_2$birth_time
  n_0 <- l_2$n_0[1]
  n_0s <- l_2$n_0

  # generic set up
  function_names <- get_function_names( # nolint internal function
    loglik_functions = loglik_functions
  )
  model_names <- get_model_names( # nolint internal function
    function_names = function_names,
    verbose = verbose
  )

  # simulate
  set.seed(seed)
  sim <- sls_sim(
    lambdas = lambdas,
    mus = mus,
    ks = ks,
    cond = cond,
    l_2 = l_2
  )
  brts <- sim$brts
  print_info(brts = brts, n_0 = n_0, cond = cond, verbose = verbose) # nolint internal function
  if (!is.list(brts)) {
    tips <- (n_0s[1] - 1) + length(brts)
  } else {
    tips <- rep(NA, length(brts))
    for (i in seq_along(brts)) {
      tips[i] <- n_0s[i] - 1 + length(brts[[i]])
    }
  }

  # maximum likelihood
  mle <- data.frame(matrix(
    NA,
    nrow = length(loglik_functions),
    # ncol must be length pars + (loglik, df, conv)
    ncol = length(sim_pars) + 3
  ))
  for (m in seq_along(loglik_functions)) {
    if (verbose == FALSE) {
      sink(tempfile())
    }
    mle_out <- sls_ml_dummy(
      loglik_function = get(function_names[m]),
      brts = brts,
      cond = cond,
      n_0 = n_0,
      start_pars = start_pars,
      optim_ids = optim_ids,
      true_pars = start_pars,
      verbose = verbose
    )
    mle[m, ] <- mle_out
    colnames(mle) <- names(mle_out)
    if (verbose == FALSE) {
      sink()
    }
  }

  # format results
  optim_ids_dummy <- c(optim_ids, FALSE, FALSE)
  results <- cbind(
    matrix(
      sim_pars,
      nrow = length(loglik_functions),
      ncol = length(sim_pars),
      byrow = TRUE
    ),
    mle,
    seed,
    cond,
    n_0,
    matrix(
      t_0s,
      nrow = length(loglik_functions),
      ncol = length(t_0s),
      byrow = TRUE
    ),
    matrix(
      tips,
      nrow = length(loglik_functions),
      ncol = length(tips),
      byrow = TRUE
    ),
    matrix(
      optim_ids_dummy,
      nrow = length(loglik_functions),
      ncol = length(optim_ids_dummy),
      byrow = TRUE
    ),
    model = model_names
  )
  if (length(t_0s) > 1) {
    t_0s_label <- paste0("t_0_", seq_len(t_0s))
  } else {
    t_0s_label <- "t_0"
  }
  if (length(tips) > 1) {
    tips_label <- paste0("tips_", seq_len(tips))
  } else {
    tips_label <- "tips"
  }
  colnames(results) <- c(
    paste0("sim_", colnames(mle[seq_len(sim_pars)])),
    colnames(mle),
    "seed",
    "cond",
    "n_0",
    t_0s_label,
    tips_label,
    paste0("optim_", names(sim_pars)),
    "model"
  )
  rownames(results) <- NULL
  results <- data.frame(results)

  # save data
  main_save_files( # nolint internal function
    project_folder = project_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed,
    sim = sim,
    results = results
  )
  print_info(brts = brts, n_0 = n_0, cond = cond, verbose = verbose) # nolint internal function
  results
}
