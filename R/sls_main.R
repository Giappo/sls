#' @title Run the full sls routine
#' @description Simulate an sls process and infer parameters maximizing the
#' likelihood(s) for given function(s)
#' @inheritParams default_params_doc
#' @details mle inference
#' @export
sls_main <- function(
  sim_pars,
  cond = 3,
  l_2 = sim_get_standard_l_2(
    crown_age = 5,
    shift_time = 2
  ),
  seed,
  start_pars = c(0.2, 0.1, 0.2, 0.1),
  optim_ids = rep(TRUE, length(start_pars)),
  loglik_functions = sls_logliks_div(),
  project_folder = NULL,
  verbose = FALSE
) {
  # check formats
  sim_pars <- as.numeric(sim_pars)
  start_pars <- as.numeric(start_pars)
  cond <- as.numeric(cond)
  seed <- as.numeric(seed)

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
    ncol = length(start_pars) + 3
  ))
  for (m in seq_along(loglik_functions)) {
    if (verbose == FALSE) {
      if (rappdirs::app_dir()$os != "win") {
        sink(tempfile())
      } else {
        sink(rappdirs::user_cache_dir())
      }
    }
    mle_out <- sls_ml(
      loglik_function = get(function_names[m]),
      brts = brts,
      cond = cond,
      n_0 = n_0,
      start_pars = start_pars,
      optim_ids = optim_ids,
      true_pars = sim_pars,
      verbose = verbose
    )
    mle[m, ] <- mle_out
    colnames(mle) <- names(mle_out)
    if (verbose == FALSE) {
      sink()
    }
  }

  # format results
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
      optim_ids,
      nrow = length(loglik_functions),
      ncol = length(optim_ids),
      byrow = TRUE
    ),
    model = model_names
  )
  if (length(t_0s) > 1) {
    t_0s_label <- paste0("t_0_", 1:length(t_0s))
  } else {
    t_0s_label <- "t_0"
  }
  if (length(tips) > 1) {
    tips_label <- paste0("tips_", 1:length(tips))
  } else {
    tips_label <- "tips"
  }
  colnames(results) <- c(
    paste0("sim_", colnames(mle[1:length(sim_pars)])),
    colnames(mle),
    "seed",
    "cond",
    "n_0",
    t_0s_label,
    tips_label,
    paste0("optim_", colnames(mle[1:length(start_pars)])),
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
