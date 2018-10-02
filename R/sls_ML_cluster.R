# rm(list = ls()); s = 1; simpars = c(0.3, 0.1, 0.6, 0.05); cond = 1; initparsopt = c(0.4, 0.15, 0.5, 0.12)

#' @title Internal sls function
#' @description Internal sls function.
#' @details This is not to be called by the user.
#' @export
sls_ML_cluster <- function(s,
                           simpars = c(0.3, 0.1, 0.6, 0.05),
                           t_d  = 4.8,
                           cond = 1,
                           initparsopt = c(0.5, 0.3, 0.5, 0.3),
                           age = 10,
                           optimmethod = 'simplex',
                           tolerance = 1E-2,
                           fun = sls::loglik_slsP)
{
  library(sls)
  fun <- eval(fun)
  s <- as.numeric(s)
  simpars <- as.numeric(simpars)
  cond <- as.numeric(cond)

  print(s)
  set.seed(s)
  # optimmethod <- 'subplex' or 'simplex'
  # pars <- c(0.3, 0.1, 0.6, 0.05)

  fun_list <- ls("package:sls")
  whichfunction1 <- NULL
  for (i in seq_along(fun_list))
  {
    if (all.equal(get(fun_list[i]), fun) == TRUE) {whichfunction1 <- i}
  }

  if (is.null(whichfunction1))
  {
    stop('This is not a likelihood function provided by sls!')
  }

  fun_name_1 <- toString(fun_list[whichfunction1]); print(fun_name_1)
  model1     <- unlist(strsplit(fun_name_1, split = 'loglik_', fixed = TRUE))[2]

  soc  <- 2
  tol  <- rep(tolerance, 3) * c(1, 10^-1, 10^-3)

  idparsopt <- c(1:2,4:5)
  idparsfix <- c(3,6,7)
  parsfix   <- c(Inf, Inf, t_d)
  idparsnoshift <- NULL
  missnumspec <- c(0,0)
  ddmodel <- 1
  maxiter <- 1000 * round((1.25)^length(idparsopt))

  pars1 <- c(simpars[1], simpars[2], Inf, simpars[3], simpars[4], Inf, t_d)
  sim   <- sls::sls_sim(pars1 = pars1, age = age, soc = soc, cond = cond)
  brtsM <- sim$brts[[1]]; brtsS <- sim$brts[[2]]; brtsM; brtsS
  NM <- (soc - 1) + length(brtsM)
  NS <- 1 + length(brtsS)
  tsplit <- min(abs(brtsM[abs(brtsM) > t_d]))
  res <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))

  pars2    <- c(res, ddmodel, cond, tsplit, 0, soc, tol, maxiter)
  outnames <- c("lambda_M", "mu_M", "K_M", "lambda_S", "mu_S", "K_S", "t_d", "LL", "df", "conv")
  simpath  <- getwd()
  datapath <- paste0(simpath, "/data")
  datafile_name <- paste0(datapath, "/sim_", s, ".RData")
  if (.Platform$OS.type != "windows" & !file.exists(datafile_name))
  {
    save(sim, file = datafile_name)
  }

  #loglik 1
  MLE_1 <- unlist(
    sls::sls_ML(loglik_function = fun,
                brtsM = brtsM,
                brtsS = brtsS,
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

  names(MLE_1) <- outnames
  out_1 <- c(MLE_1, NM, NS, s);
  names(out_1) <- c(names(MLE_1), "tips_M", "tips_S", "tree_id")
  file_name <- paste0(simpath, "/", model1, "_MLE", s, ".txt")
  if (!file.exists(file_name))
  {
    write.table(matrix(out_1, ncol = length(out_1)), file = file_name,
                append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  }

  return(out_1 = out_1)
}

#' @title Internal sls function
#' @description Internal sls function.
#' @details This is not to be called by the user.
#' @export
sls_ML_cluster2 <- function(s,
                           simpars = c(0.3, 0.1, 0.6, 0.05),
                           t_d  = 4.8,
                           cond = 1,
                           initparsopt = c(0.5, 0.3, 0.5, 0.3),
                           age = 10,
                           optimmethod = 'simplex',
                           tolerance = 1E-2,
                           fun1 = sls::loglik_slsP,
                           fun2 = sls::loglik_DDD)
{
  set.seed(s)
  # optimmethod <- 'subplex' or 'simplex'
  # pars <- c(0.3, 0.1, 0.6, 0.05)

  fun_list <- ls("package:sls")
  for (i in seq_along(fun_list))
  {
    if (all.equal(get(fun_list[i]), fun1) == TRUE) {whichfunction1 <- i}
    if (all.equal(get(fun_list[i]), fun2) == TRUE) {whichfunction2 <- i}
  }

  fun_name_1 <- toString(fun_list[whichfunction1])
  fun_name_2 <- toString(fun_list[whichfunction2])
  model1 <- unlist(strsplit(fun_name_1, split = 'loglik_', fixed = TRUE))[2]
  model2 <- unlist(strsplit(fun_name_2, split = 'loglik_', fixed = TRUE))[2]

  soc  <- 2
  tol  <- rep(tolerance, 3) * c(1, 10^-1, 10^-3)

  idparsopt <- c(1:2,4:5)
  idparsfix <- c(3,6,7)
  parsfix <- c(Inf, Inf, t_d)
  idparsnoshift <- NULL
  missnumspec <- c(0,0)
  ddmodel <- 1
  maxiter <- 1000 * round((1.25)^length(idparsopt))

  pars1 <- c(simpars[1], simpars[2], Inf, simpars[3], simpars[4], Inf, t_d)
  sim   <- sls::sls_sim(pars1 = pars1, age = age, soc = soc, cond = cond)
  brtsM <- sim$brts[[1]]; brtsS <- sim$brts[[2]]; brtsM; brtsS
  NM <- (soc - 1) + length(brtsM)
  NS <- 1 + length(brtsS)
  tsplit <- min(abs(brtsM[abs(brtsM) > t_d]))
  res <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))

  pars2    <- c(res, ddmodel, cond, tsplit, 0, soc, tol, maxiter)
  outnames <- c("lambda_M", "mu_M", "K_M", "lambda_S", "mu_S", "K_S", "t_d", "LL", "df", "conv")
  simpath  <- getwd()
  datapath <- paste0(simpath, "/data")
  if (.Platform$OS.type != "windows")
  {
    save(sim, file = paste0(datapath, "/sim_", s, ".RData"))
  }

  #loglik 1
  MLE_1 <- unlist(
    sls::sls_ML(loglik_function = fun1,
                brtsM = brtsM,
                brtsS = brtsS,
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

  names(MLE_1) <- outnames
  out_1 <- c(MLE_1, NM, NS, s);
  names(out_1) <- c(names(MLE_1), "tips_M", "tips_S", "tree_id")
  write.table(matrix(out_1, ncol = length(out_1)), file = paste0(simpath, "/", model1,"_MLE", s, ".txt"),
              append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")

  #DDD_KI
  MLE_2 <- unlist(
    sls::sls_ML(loglik_function = fun2,
                brtsM = brtsM,
                brtsS = brtsS,
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
  out_2 <- c(MLE_2, NM, NS, s);
  names(out_2) <- c(names(MLE_2), "tips_M", "tips_S", "tree_id")
  write.table(matrix(out_2, ncol = length(out_2)), file = paste0(simpath, "/", model2,"_MLE", s, ".txt"),
              append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  return(list(out_1 = out_1, out_2 = out_2))
}
