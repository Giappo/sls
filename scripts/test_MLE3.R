rm(list = ls()); s = 4
simpars <- c(0.5, 0.2, 0.6, 0.1) #c(0.4, 0.2, 0.6, 0.1)
initparsopt <- c(0.5, 0.3, 0.5, 0.3)
cond <- 1
t_d <- 4
tolerance <- 1e-2
soc <- 2
missnumspec <- c(0,0)
res <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
# t_sub <- system.time(test_sub <- sls::sls_ML_cluster(s = s,
#                             simpars = simpars,
#                             cond = cond,
#                             initparsopt = c(0.5, 0.3, 0.5, 0.3),
#                             t_d = t_d,
#                             optimmethod = 'subplex',
#                             tolerance = 1e-1)
# )
set.seed(s)
pars1 <- c(simpars[1], simpars[2], Inf, simpars[3], simpars[4], Inf, t_d)
brts <- sls_sim2(pars1 = pars1, age = 10, soc = soc, cond = cond)$brts; brtsM <- brts[[1]]; brtsS <- brts[[2]]; tsplit <- max(brtsM[brtsM < min(brtsS)])
t_lik <- system.time(test_lik <- sls::lik_shift_P2(pars1 = pars1, pars2 = c(res, 1, cond, tsplit, 0 , soc),
                                                   brtsM = brtsM, brtsS = brtsS))
t_sim <- system.time(test_sim <- sls::sls_ML_cluster(s = s,
                            simpars = simpars,
                            cond = cond,
                            initparsopt = initparsopt,
                            t_d = t_d,
                            optimmethod = 'simplex',
                            tolerance = 1e-2)
)
print(test_sim); print(simpars)
# print(test_sub); print(test_sim);
# MLE <- list(maxs <- 10)
# for(s in 1:maxs) {MLE[[s]] <- unlist(sls::sls_ML_cluster(s = s))}
