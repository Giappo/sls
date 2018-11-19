rm(list = ls()); library(DDD)

#1 #typo for pars1[4]
?dd_KI_loglik

#2 # in dd_KI_sim subclade exhibits species born before the time of the shift
set.seed(3)
pars <- c(0.3, 0.1, 0.6, 0.05) #lambdaM, muM, lambdaS, muS
t0 <- c(10, 4) #starting time for each clade
test <- DDD::dd_KI_sim(pars = c(pars[1], pars[2], Inf, pars[3], pars[4], Inf, t0[2]), age = t0[1], ddmodel = 1)
L <- test$L
L0 <- L[L[,5] == 0, 1:5]; dim(L0) <- c(sum(L[,5] == 0), 5); head(L0) #selecting only the main clade L[,5]==0
DDD::L2brts(L0) #why it doesn't work?
head(L) #how come that one of the two crown species belongs to the subclade if the shift is at t0[2] = 4?
L1 <- L[L[,5] == 1, 1:5]; dim(L1) <- c(sum(L[,5] == 1), 5); head(L1) #selecting only the sub clade L[,5]==0

#3 can L2brts work with a phylogeny starting with a stem?
l1 <- c(10, 0, -1, -1)
l2 <- c(8, -1, -2, -1)
l3 <- c(5, -2, -3, -1)
L  <- rbind(l1, l2, l3)
DDD::L2brts(L)

#4
pars <- c(0.3, 0.1, 0.6, 0.05) #lambdaM, muM, lambdaS, muS
DDD::dd_KI_loglik(pars1 = c(pars[1], pars[2], Inf, pars[3], pars[4], Inf, 5),
                  pars2 = c(100, 1, 1, 5, 0, 2),
                  brtsM = brtsM,
                  brtsS = NULL,
                  missnumspec = c(0, 0))


