# rm(list = ls()); seed = 1; sim_pars = c(0.3, 0.1, 0.6, 0.05); cond = 1; initparsopt = c(0.4, 0.15, 0.5, 0.12) # nolint test parameters

#' @title Internal sls function
#' @description Internal sls function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
sls_ml_cluster <- function(
  seed,
  sim_pars = c(0.3, 0.1, 0.6, 0.05),
  t_d  = 4.8,
  cond = 1,
  initparsopt = c(0.5, 0.3, 0.5, 0.3),
  age = 10,
  optimmethod = "simplex",
  tolerance = 1E-2,
  fun = sls::loglik_sls_p
) {
  library(sls)
  fun <- eval(fun)
  seed <- as.numeric(seed)
  sim_pars <- as.numeric(sim_pars)
  cond <- as.numeric(cond)

  print(seed)
  set.seed(seed)
  # optimmethod <- 'subplex' or 'simplex' # nolint
  # pars <- c(0.3, 0.1, 0.6, 0.05) # nolint

  fun_list <- ls("package:sls")
  whichfunction1 <- NULL
  for (i in seq_along(fun_list)) {
    if (all.equal(get(fun_list[i]), fun) == TRUE) {
      whichfunction1 <- i
    }
  }

  if (is.null(whichfunction1)) {
    stop("This is not a likelihood function provided by sls!")
  }

  fun_name_1 <- toString(fun_list[whichfunction1]); print(fun_name_1)
  model1     <- unlist(strsplit(fun_name_1, split = "loglik_", fixed = TRUE))[2]

  soc  <- 2
  tol  <- rep(tolerance, 3) * c(1, 1e-1, 1e-3)

  idparsopt <- c(1:2, 4:5)
  idparsfix <- c(3, 6, 7)
  parsfix   <- c(Inf, Inf, t_d)
  idparsnoshift <- NULL

  pars1 <- c(sim_pars[1], sim_pars[2], Inf, sim_pars[3], sim_pars[4], Inf, t_d)
  sim   <- sls::sls_sim(pars1 = pars1, age = age, soc = soc, cond = cond)
  brts_m <- sim$brts[[1]]; brts_s <- sim$brts[[2]]; brts_m; brts_s
  NM <- (soc - 1) + length(brts_m)
  NS <- 1 + length(brts_s)
  tsplit <- min(abs(brts_m[abs(brts_m) > t_d]))

  # ddmodel <- 1 # nolint
  # maxiter <- 1000 * round(1.25 ^ length(idparsopt)) # nolint
  # res <- 10 * (1 + length(c(brts_m, brts_s)) + sum(missnumspec)) # nolint
  # pars2    <- c(res, ddmodel, cond, tsplit, 0, soc, tol, maxiter) # nolint
  outnames <- c(
    "lambda_M",
    "mu_M",
    "K_M",
    "lambda_S",
    "mu_S",
    "K_S",
    "t_d",
    "LL",
    "df",
    "conv"
  )
  simpath  <- getwd()
  datapath <- file.path(simpath, "data")
  datafile_name <- file.path(
    datapath,
    paste0("sim_", seed, ".RData")
  )
  if (.Platform$OS.type != "windows" & !file.exists(datafile_name)) {
    if (!file.exists(datapath)) {
      dir.create(datapath, showWarnings = FALSE)
    }
    save(sim, file = datafile_name)
  }

  file_name <- paste0(simpath, "/", model1, "_MLE", seed, ".txt")
  if (!file.exists(file_name)) {
    MLE_1 <- unlist(
      sls::sls_ml(
        loglik_function = fun,
        brts_m = brts_m,
        brts_s = brts_s,
        tsplit = tsplit,
        idparsopt = idparsopt,
        initparsopt = initparsopt,
        idparsfix = idparsfix,
        parsfix = parsfix,
        idparsnoshift = idparsnoshift,
        cond = cond,
        tol = tol,
        optimmethod = optimmethod
      )
    )

    names(MLE_1) <- outnames
    out_1 <- c(MLE_1, NM, NS, seed);
    names(out_1) <- c(names(MLE_1), "tips_M", "tips_S", "tree_id")

    utils::write.table(matrix(out_1, ncol = length(out_1)), file = file_name,
                append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  }

  return(out_1 = out_1)
}

#' @title Internal sls function
#' @description Internal sls function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
sls_ml_cluster2 <- function(
  seed,
  sim_pars = c(0.3, 0.1, 0.6, 0.05),
  t_d  = 4.8,
  cond = 1,
  initparsopt = c(0.5, 0.3, 0.5, 0.3),
  age = 10,
  optimmethod = "simplex",
  tolerance = 1E-2,
  fun1 = sls::loglik_sls_p,
  fun2 = sls::loglik_ddd
) {
  set.seed(seed)
  # optimmethod <- 'subplex' or 'simplex' #nolint
  # pars <- c(0.3, 0.1, 0.6, 0.05) # nolint

  fun_list <- ls("package:sls")
  for (i in seq_along(fun_list)) {
    if (all.equal(get(fun_list[i]), fun1) == TRUE) {
      whichfunction1 <- i
    }
    if (all.equal(get(fun_list[i]), fun2) == TRUE) {
      whichfunction2 <- i
    }
  }

  fun_name_1 <- toString(fun_list[whichfunction1])
  fun_name_2 <- toString(fun_list[whichfunction2])
  model1 <- unlist(strsplit(fun_name_1, split = "loglik_", fixed = TRUE))[2]
  model2 <- unlist(strsplit(fun_name_2, split = "loglik_", fixed = TRUE))[2]

  soc  <- 2
  tol  <- rep(tolerance, 3) * c(1, 1e-1, 1e-3)

  idparsopt <- c(1:2, 4:5)
  idparsfix <- c(3, 6, 7)
  parsfix <- c(Inf, Inf, t_d)
  idparsnoshift <- NULL

  pars1 <- c(sim_pars[1], sim_pars[2], Inf, sim_pars[3], sim_pars[4], Inf, t_d)
  sim   <- sls::sls_sim(pars1 = pars1, age = age, soc = soc, cond = cond)
  brts_m <- sim$brts[[1]]; brts_s <- sim$brts[[2]]; brts_m; brts_s
  NM <- (soc - 1) + length(brts_m)
  NS <- 1 + length(brts_s)
  tsplit <- min(abs(brts_m[abs(brts_m) > t_d]))

  # ddmodel <- 1 # nolint
  # maxiter <- 1000 * round(1.25 ^ length(idparsopt)) # nolint
  # res <- 10 * (1 + length(c(brts_m, brts_s)) + sum(missnumspec)) # nolint
  # pars2    <- c(res, ddmodel, cond, tsplit, 0, soc, tol, maxiter) # nolint
  outnames <- c(
    "lambda_M",
    "mu_M",
    "K_M",
    "lambda_S",
    "mu_S",
    "K_S",
    "t_d",
    "LL",
    "df",
    "conv"
  )
  simpath  <- getwd()
  datapath <- file.path(simpath, "data")
  if (.Platform$OS.type != "windows") {
    filename <- file.path(
      datapath,
      paste0("sim_", seed, ".RData")
    )
    save(sim, file = filename)
  }

  #loglik 1
  MLE_1 <- unlist(
    sls::sls_ml(
      loglik_function = fun1,
      brts_m = brts_m,
      brts_s = brts_s,
      tsplit = tsplit,
      idparsopt = idparsopt,
      initparsopt = initparsopt,
      idparsfix = idparsfix,
      parsfix = parsfix,
      idparsnoshift = idparsnoshift,
      cond = cond,
      tol = tol,
      optimmethod = optimmethod
    )
  )

  names(MLE_1) <- outnames
  out_1 <- c(MLE_1, NM, NS, seed);
  names(out_1) <- c(names(MLE_1), "tips_M", "tips_S", "tree_id")
  utils::write.table(
    matrix(
      out_1,
      ncol = length(out_1)
    ),
    file = paste0(simpath, "/", model1, "_MLE", seed, ".txt"),
    append = TRUE,
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
  )

  #DDD_KI
  MLE_2 <- unlist(
    sls::sls_ml(
      loglik_function = fun2,
      brts_m = brts_m,
      brts_s = brts_s,
      tsplit = tsplit,
      idparsopt = idparsopt,
      initparsopt = initparsopt,
      idparsfix = idparsfix,
      parsfix = parsfix,
      idparsnoshift = idparsnoshift,
      cond = cond,
      tol = tol,
      optimmethod = optimmethod)
  )

  names(MLE_2) <- outnames
  out_2 <- c(MLE_2, NM, NS, seed);
  names(out_2) <- c(names(MLE_2), "tips_M", "tips_S", "tree_id")
  utils::write.table(
    matrix(out_2, ncol = length(out_2)),
    file = paste0(simpath, "/", model2, "_MLE", seed, ".txt"),
    append = TRUE,
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
  )
  return(list(
    out_1 = out_1,
    out_2 = out_2)
  )
}
