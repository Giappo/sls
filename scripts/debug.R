rm(list = ls())
#data
brtsM      <- c(-10, -10)
brtsS      <- c(-3)
lambdas    <- c(0.5, 0.3)
mus        <- c(0.4, 0.1)

#TEST LIK
for (cond in 0:4) {
  testP <- sls::lik_shift_P(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus, cond = cond)
  testQ <- sls::lik_shift_Q(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus, cond = cond)
  testit::assert(abs(1 - testP/testQ) <= 1e-5)
}

#TEST Pc1
#SLS_Pc1 = PM1 * (PS * PM2cs + (1 - PS) * PM2cp)
#PM1 = survival of one branch from crown to present
#PS = survival of the subclade from shift to present
#PM2cs = survival of second branch from crown to shift
#PM2cp = survival of second branch from crown to present
# SLS <- sls::debug_Pc_1shift(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus); SLS
SLS <- sls::Pc_1shift(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus); SLS
#DDD_Pc1 = 2 * (PS * PM2 + PM12)
#PS = survival of the subclade from shift to present
#PM2 = survival of second branch from crown to present
#PM12 = survival of both crown species to present
DDD <- sls::debug_DDD_conditioning(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus, lx = 200); DDD
# DDD_Pc1 <- sls::DDD_conditioning(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus, lx = 200); DDD

SLS$Pc1
DDD_Pc1 <- 2 * (DDD$PS * DDD$PM2 + DDD$PM12); DDD_Pc1
SLS$Pc2
DDD_Pc2 <- 2 * DDD$PS * DDD$PM12; DDD_Pc2
SLS$Pc3
DDD_Pc3 <- 2 * (1 - DDD$PS) * DDD$PM12; DDD_Pc3

testit::assert(abs(SLS_Pc1 - DDD_Pc1) < 1e-10) #the two Pc1 don't match

#there is a strong possibility that comparing DDD_Pc1 with a term that comprises (1-PS) is not a right comparison
SLS_Pc1_alt <- SLS$PM1 * (SLS$PS * SLS$PM2cs + SLS$PM2cp); SLS_Pc1_alt
testit::assert(abs(SLS_Pc1_alt - DDD_Pc1) < 1e-10) #this doesn't work either

#PS is the same for both approaches. The error is in the main clade
testit::assert(abs(DDD$PS - SLS$PS) < 1e-10) #survival of the subclade: they match

testit::assert(abs(DDD$PM2cs - SLS$PM2cs) < 1e-10)
testit::assert(abs(SLS$PM1 * SLS$PM2cp - 2 * DDD$PM12) < 1e-10) #survival of both mainclade's crown descendants until the present. They don't match

#TEST Pc3
#Pc3 should make both main branch survive as well as the subclade
SLS_Pc3 <- SLS$PM1 * SLS$PS * SLS$PM2cs; SLS_Pc3
DDD_Pc3 <- 2 * DDD$PS * DDD$PM12; DDD_Pc3
testit::assert(abs(SLS_Pc3 - DDD_Pc3) < 1e-10) #they don't match :(
