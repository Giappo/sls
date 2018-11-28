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

  model_names <- which_function <- rep(NA, length(models))
  for (m in seq_along(models)) {
    fun <- eval(models[m])[[1]]
    fun_list <- ls("package:sls")
    if (is.character(models[m])) {
      which_function[m] <- which(fun_list == models[m])
    } else {
      for (i in seq_along(fun_list)) {

        if (all.equal(get(fun_list[i]), fun) == TRUE) {
          which_function[m] <- i
        }
      }
    }
    if (is.null(which_function[m])) {
      stop("This is not a likelihood function provided by sls!")
    }
    fun_name_1 <- toString(fun_list[which_function[m]])
    model_names[m] <- unlist(strsplit(
      fun_name_1,
      split = "loglik_",
      fixed = TRUE
    ))[2]
  }
  if (verbose == TRUE) {
   cat("You are using the functions:", model_names)
  }

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
    mle <- sls_ml(
      loglik_function = get(fun_list[which_function[m]]),
      brts_m = sim$brts[[1]],
      brts_s = sim$brts[[2]],
      start_pars = start_pars,
      cond = cond,
      n_0 = l_2$n_0[1],
      verbose = FALSE
    )
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
    sim_path  <- system.file("extdata", package = "sls")
    if (!file.exists(sim_path)) {
      dir.create(sim_path, showWarnings = FALSE)
    }
  } else {
    sim_path  <- getwd()
  }
  data_path <- file.path(sim_path, "data")
  data_file_name <- file.path(
    data_path,
    paste0("sim_", seed, ".RData")
  )
  if (!file.exists(data_file_name)) {
    if (!file.exists(data_path)) {
      dir.create(data_path, showWarnings = FALSE)
    }
  }
  save(sim, file = data_file_name)
  results_file_name <- file.path(sim_path, paste0("sls_mle", seed, ".txt"))
  utils::write.csv(
    x = out,
    file = results_file_name
  )

  out
}
