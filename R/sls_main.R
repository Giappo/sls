#' @title Run the full sls routine
#' @description Simulate an sls process and infer parameters maximizing the
#' likelihood(s) for given function(s)
#' @inheritParams default_params_doc
#' @details mle inference
#' @export
sls_main <- function(
  seed,
  sim_pars,
  cond,
  l_2 = sim_get_standard_l_2(
    crown_age = 5,
    shift_time = 2
  ),
  start_pars = c(0.2, 0.1, 0.2, 0.1),
  models = sls_logliks_div(),
  verbose = TRUE
) {
  # set up
  lambdas <- sim_pars[c(1, 3)]
  mus <- sim_pars[c(2, 4)]
  ks <- c(Inf, Inf)
  pkg_name <- sls_pkg_name()

  function_names <- sls_get_function_names(
    models = models
  )
  model_names <- sls_get_model_names(
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
  tips_m <- (l_2[l_2$clade_id == 1, ]$n_0 - 1) + length(sim$brts[[1]])
  tips_s <- (l_2[l_2$clade_id == 2, ]$n_0 - 1) + length(sim$brts[[2]])

  # maximum likelihood
  results <- data.frame(matrix(
    NA,
    nrow = length(models),
    # ncol must be length pars + (loglik, df, conv) + (tips_m, tips_s, seed)
    ncol = length(start_pars) + 3 + 3
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
      brts_m = sim$brts[[1]],
      brts_s = sim$brts[[2]],
      start_pars = start_pars,
      cond = cond,
      n_0 = l_2$n_0[1],
      verbose = FALSE
    )
    if (verbose == FALSE) {
      sink()
    }
    results[m, ] <- data.frame(
      mle,
      tips_m = tips_m,
      tips_s = tips_s,
      seed = seed
    )
  }

  # format output
  colnames(results) <- c(
    colnames(mle),
    "tips_m",
    "tips_s",
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
  out <- data.frame(out)

  # save data
  if (.Platform$OS.type == "windows") {
    sim_path  <- system.file("extdata", package = pkg_name)
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

  out
}
