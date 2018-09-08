# rm(list = ls()); s = 1; simpars = c(0.3, 0.1, 0.6, 0.05); cond = 1; initparsopt = c(0.4, 0.15, 0.5, 0.12)

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
#' @export
sls_ML_cluster <- function(s,
                           simpars = c(0.3, 0.1, 0.6, 0.05),
                           t_d  = 4.8,
                           cond = 1,
                           initparsopt = c(0.5, 0.3, 0.5, 0.3),
                           age = 10,
                           optimmethod = 'simplex',
                           tolerance = 1E-2)
{
  # optimmethod <- 'subplex' #'simplex'
  set.seed(s)
  # pars <- c(0.3, 0.1, 0.6, 0.05)

  soc  <- 2
  tol  <- rep(tolerance, 3) * c(1, 10^-1, 10^-3)

  idparsopt <- c(1:2,4:5)
  idparsfix <- c(3,6,7)
  parsfix <- c(Inf, Inf, t_d)
  idparsnoshift <- NULL
  missnumspec <- c(0,0)
  ddmodel <- 1
  maxiter <- 1000 * round((1.25)^length(idparsopt))

  # sim <- sls::sls_sim(pars = simpars, t0 = c(age, t_d), starting_species = c(soc, 1), cond = cond)
  pars1 <- c(simpars[1], simpars[2], Inf, simpars[3], simpars[4], Inf, t_d)
  sim   <- sls::sls_sim2(pars1 = pars1, age = age, soc = soc, cond = cond)
  brtsM <- sim$brts[[1]]; brtsS <- sim$brts[[2]]; brtsM; brtsS
  # tsplit <- abs(max(brtsM[abs(brtsM) > t_d]))
  tsplit <- min(abs(brtsM[abs(brtsM) > t_d]))
  res <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))

  # pars1a <- rep(NA, 7); pars1a[idparsopt] <- initparsopt; pars1a[idparsfix] <- parsfix
  # sim2 <- DDD::dd_KI_sim(pars = pars1a, age = age, ddmodel = ddmodel)

  # pars1 <- c(pars[1], pars[2], Inf, pars[3], pars[4], Inf, t_d)
  pars2 <- c(res, ddmodel, cond, tsplit, 0, soc, tol, maxiter)
  outnames <- c("lambda_M", "mu_M", "K_M", "lambda_S", "mu_S", "K_S", "t_d", "LL", "df", "conv")
  simpath  <- getwd()
  datapath <- paste0(simpath, "/data")
  save(sim, file = paste0(datapath, "/sim_", s, ".RData"))

  #sls
  MLE_sls <- unlist(
    sls::sls_ML2(loglik_function = sls:::lik_shift_P2,
                 brtsM = brtsM, brtsS = brtsS, tsplit = tsplit, idparsopt = idparsopt,
                 initparsopt = initparsopt,
                 idparsfix = idparsfix,
                 parsfix = parsfix,
                 idparsnoshift = idparsnoshift,
                 cond = cond,
                 tol = tol,
                 optimmethod = optimmethod)
  )

  names(MLE_sls) <- outnames
  out_sls  <- c(MLE_sls, s);
  names(out_sls) <- c(names(MLE_sls), "tree id")
  # print(out)
  # sink()
  # print(out)
  write.table(matrix(out_sls, ncol = length(out_sls)), file = paste0(simpath, "/sls_MLE", s, ".txt"),
              append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")

  #DDD_KI
  MLE_DDD <- unlist(
    sls::sls_ML2(loglik_function = sls:::lik_shift_DDD2,
                 brtsM = brtsM, brtsS = brtsS, tsplit = tsplit, idparsopt = idparsopt,
                 initparsopt = initparsopt,
                 idparsfix = idparsfix,
                 parsfix = parsfix,
                 idparsnoshift = idparsnoshift,
                 cond = cond,
                 tol = tol,
                 optimmethod = optimmethod)
  )
  names(MLE_DDD) <- outnames
  out_DDD  <- c(MLE_DDD, s);
  names(out_DDD) <- c(names(MLE_DDD), "tree id")
  # print(out)
  # sink()
  # print(out)
  write.table(matrix(out_DDD, ncol = length(out_DDD)), file = paste0(simpath, "/DDD_MLE", s, ".txt"),
              append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  return(list(out_sls = out_sls, out_DDD = out_DDD))
}
