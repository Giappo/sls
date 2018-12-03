#' @title Run the full sls routine
#' @description Simulate an sls process and infer parameters maximizing the
#' likelihood(s) for given function(s)
#' @inheritParams default_params_doc
#' @details mle inference
#' @export
sls_main <- function(
  seed,
  sim_pars,
  cond = 3,
  l_2 = sim_get_standard_l_2(
    crown_age = 5,
    shift_time = 2
  ),
  start_pars = c(0.2, 0.1, 0.2, 0.1),
  models = sls_logliks_div(),
  verbose = FALSE
) {
  # specific set up
  lambdas <- sim_pars[c(1, 3)]
  mus <- sim_pars[c(2, 4)]
  ks <- c(Inf, Inf)
  n_0s <- l_2$n_0

  # generic set up
  pkg_name <- get_pkg_name() # nolint internal function
  function_names <- get_function_names( # nolint internal function
    models = models
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
  print_info(brts = brts, n_0s = n_0s, cond = cond, verbose = verbose) # nolint internal function
  if (!is.list(brts)) {
    tips <- (n_0s[1] - 1) + length(brts)
  } else {
    tips <- rep(NA, length(brts))
    for (i in seq_along(brts)) {
      tips[i] <- n_0s[i] - 1 + length(brts[[i]])
    }
  }

  # maximum likelihood
  results <- data.frame(matrix(
    NA,
    nrow = length(models),
    # ncol must be length pars + (loglik, df, conv, seed) + length(tips)
    ncol = length(start_pars) + 4 + length(tips)
  ))
  for (m in seq_along(models)) {
    if (verbose == FALSE) {
      if (rappdirs::app_dir()$os != "win") {
        sink(file.path(rappdirs::user_cache_dir(), "ddd"))
      } else {
        sink(rappdirs::user_cache_dir())
      }
    }
    mle <- sls_ml(
      loglik_function = get(function_names[m]),
      brts = sim$brts,
      start_pars = start_pars,
      cond = cond,
      n_0 = l_2$n_0[1],
      verbose = verbose
    )
    if (verbose == FALSE) {
      sink()
    }
    dim(tips) <- c(1, length(tips))
    results[m, ] <- data.frame(
      cbind(mle, tips, seed)
    )
  }

  # format output
  colnames(results) <- c(
    colnames(mle),
    paste0("tips_", 1:length(tips)),
    "seed"
  )
  out <- cbind(
    matrix(
      sim_pars,
      nrow = length(models),
      ncol = length(start_pars),
      byrow = TRUE
    ),
    results,
    model = model_names
  )
  colnames(out) <- c(
    paste0("sim_", colnames(mle[1:length(start_pars)])),
    colnames(results),
    "model"
  )
  rownames(out) <- NULL
  out <- data.frame(out)

  # save data
  if (.Platform$OS.type == "windows") {
    if (!("extdata" %in% list.files(system.file(package = pkg_name)))) {
      dir.create(file.path(
        system.file(package = pkg_name),
        "extdata"
      ))
    }
    sim_path <- system.file("extdata", package = pkg_name)
    if (!file.exists(sim_path)) {
      dir.create(sim_path, showWarnings = FALSE)
    }
  } else {
    sim_path  <- getwd()
  }
  data_path <- file.path(sim_path, "data")
  data_file_name <- file.path(
    data_path,
    paste0(pkg_name, "_sim_", seed, ".RData")
  )
  if (!file.exists(data_file_name)) {
    if (!file.exists(data_path)) {
      dir.create(data_path, showWarnings = FALSE)
    }
  }
  save(sim, file = data_file_name)
  results_file_name <- file.path(
    sim_path,
    paste0(pkg_name, "_mle_", seed, ".txt")
  )
  utils::write.csv(
    x = out,
    file = results_file_name
  )
  print_info(brts = brts, n_0s = n_0s, cond = cond, verbose = verbose) # nolint internal function
  out
}
