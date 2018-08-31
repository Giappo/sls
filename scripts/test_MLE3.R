rm(list = ls()); s = 4
simpars <- c(0.5, 0.2, 0.6, 0.1) #c(0.4, 0.2, 0.6, 0.1)
initparsopt <- c(0.5, 0.3, 0.5, 0.3)
cond <- 1
age <- 10
t_d <- 4
tolerance <- 1e-2
soc <- 2
missnumspec <- c(0,0)
set.seed(s)
pars1 <- c(simpars[1], simpars[2], Inf, simpars[3], simpars[4], Inf, t_d)
brts  <- sls::sls_sim2(pars1 = pars1, age = age, soc = soc, cond = cond)$brts; brtsM <- brts[[1]]; brtsS <- brts[[2]]; tsplit <- max(brtsM[brtsM < min(brtsS)])
res   <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
pars2 <- c(res, 1, cond, tsplit, 0 , soc)
t_lik <- system.time(test_lik_sls <- sls::lik_shift_P2(pars1 = pars1, pars2 = pars2,
                                                       brtsM = brtsM, brtsS = brtsS)); test_lik_sls
t_lik <- system.time(test_lik_DDD <- sls::lik_shift_DDD2(pars1 = pars1, pars2 = pars2,
                                                         brtsM = brtsM, brtsS = brtsS)); test_lik_DDD
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
