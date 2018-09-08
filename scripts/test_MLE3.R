rm(list = ls()); s = 10 #10 and 11 give brtsS == NULL and DDD fails on that

for (s in 1:9){
simpars <- c(0.3, 0.1, 0.6, 0.08) #c(0.4, 0.2, 0.6, 0.1)
initparsopt <- c(0.5, 0.3, 0.5, 0.3)
cond <- 1
age <- 10
tolerance <- 1e-2
soc <- 2
missnumspec <- c(0,0)
t_d <- 4.8
set.seed(s)
pars1 <- c(initparsopt[1], initparsopt[2], Inf, initparsopt[3], initparsopt[4], Inf, t_d)
sim   <- sls::sls_sim2(pars1 = pars1, age = age, soc = soc, cond = cond)
brts  <- sim$brts; brtsM <- brts[[1]]; brtsS <- brts[[2]]
tsplit <- abs(max(brtsM[abs(brtsM) > t_d]))
res   <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
pars2 <- c(res, 1, cond, tsplit, 0 , soc)

# t_sls <- system.time(test_lik_sls <- sls::lik_shift_P2(pars1 = pars1, pars2 = pars2,
#                                                        brtsM = brtsM, brtsS = brtsS)); test_lik_sls
# t_DDD <- system.time(test_lik_DDD <- sls::lik_shift_DDD2(pars1 = pars1, pars2 = pars2,
#                                                          brtsM = brtsM, brtsS = brtsS)); test_lik_DDD
t_MLE <- system.time(test_MLE <- sls::sls_ML_cluster(s = s,
                                                     simpars = simpars,
                                                     cond = cond,
                                                     initparsopt = initparsopt,
                                                     t_d = t_d,
                                                     optimmethod = 'simplex',
                                                     tolerance = 1e-2)
)
print(test_MLE); print(simpars)
}
# print(test_sub); print(test_sim);
# MLE <- list(maxs <- 10)
# for(s in 1:maxs) {MLE[[s]] <- unlist(sls::sls_ML_cluster(s = s))}
# load(paste0("F:\Dropbox\University\Progress\RQ4-single_lineage_rate_shifts\sls\data\sim_1.Rdata"))
# load(paste0("F://Dropbox//University//Progress//RQ4-single_lineage_rate_shifts//sls//data//sim_1.Rdata"))


# DDD::dd_KI_loglik(pars1 = pars1, pars2 = pars2,
                  # brtsM = brtsM, brtsS = brtsS, missnumspec = missnumspec)
#
#
# DDD::dd_KI_ML(brtsM = brtsM, brtsS = brtsS, tsplit = tsplit, initparsopt = initparsopt,
#               parsfix = parsfix, idparsopt = idparsopt, idparsfix = idparsfix, idparsnoshift = NULL, res = res,
#               cond = cond, soc = soc)
