rm(list = ls())

pars <- c(0.3, 0.1, 0.6, 0.05)
t_d <- 4.8
cond <- 1
soc <- 2
idparsopt = c(1:2,4:5)
initparsopt = c(0.6, 0.2, 0.3, 0.2)
idparsfix = c(3,6,7)
parsfix = c(Inf, Inf, t_d)
idparsnoshift = NULL
missnumspec <- c(0,0)
ddmodel <- 1
tol = c(1E-3, 1E-3, 1E-3)
maxiter = 1000 * round((1.25)^length(idparsopt))
pars.transform <- 1

maxs <- 10
MLE <- matrix(NA, nrow = maxs, ncol = 10); s <- 1
for (s in 1:maxs)
{
  set.seed(s)
  sim <- sls::sls_sim(pars = pars, t0 = c(10, t_d), starting_species = c(soc, 1), cond = cond)
  brtsM <- sim$brts[[1]]; brtsS <- sim$brts[[2]]; brtsM; brtsS
  res = 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
  tsplit <- max(brtsM[brtsM < min(brtsS)])
  pars1 <- c(pars[1], pars[2], Inf, pars[3], pars[4], Inf, t_d)
  pars2 = c(res, ddmodel, cond, tsplit, 0, soc, tol, maxiter)
  MLE[s,] <- unlist(
    sls::sls_ML2(brtsM = brtsM, brtsS = brtsS, tsplit = tsplit, idparsopt = idparsopt,
                           initparsopt = initparsopt,
                           idparsfix = idparsfix,
                           parsfix = parsfix,
                           idparsnoshift = idparsnoshift,
                           cond = cond,
                           tol = tol,
                           optimmethod = 'simplex')
  )
}; colnames(MLE) <- c("lambda_M", "mu_M", "K_M", "lambda_S", "mu_S", "K_S", "t_d", "LL", "df", "conv"); MLE
MLE2 <- MLE[MLE[,1] != -1 & !is.na(MLE[,1]),]; head(MLE2)
MLE_medians <- apply(MLE2, MARGIN = 2, FUN = function(x) median(x))[1:7]; MLE_medians
pars1



# sls::lik_shift_P2(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS, missnumspec = missnumspec)

# DDD:::dd_KI_ML(brtsM = brtsM, brtsS = brtsS, tsplit = tsplit, idparsopt = idparsopt,
#                initparsopt = initparsopt,
#                idparsfix = idparsfix,
#                parsfix = parsfix,
#                idparsnoshift = idparsnoshift,
#                cond = cond,
#                tol = tol,
#                optimmethod = 'simplex')
