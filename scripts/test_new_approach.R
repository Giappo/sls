rm(list = ls())
# #conditional probabilities
# #Pc1 = AMNAT conditioning: PM1 x (PS + (1 - PS) * PM2)
# #Pc2 = alternative: PM1 x PM2 x PS
# #PC3 = "real": PM1 x PM2(til the shift) x PS
# brtsM      <- c(-10, -10)
# brtsS      <- c(-4)
# lambdas    <- c(0.5, 0.4)
# mus        <- c(0.4, 0.2)
# Pc1 <- sls:::Pc1_1shift(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus); Pc1
# Pc3 <- sls:::Pc3_1shift(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus); Pc3
# Qc3 <- sls:::Qc3_1shift(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus, maxN = 50); Qc3
# Qc30 <- sls:::Qc3_1shift0(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus, lx = 300); Qc30
# Pc1_DDD <- sls:::DDD_conditioning(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus, lx = 100); Pc1_DDD

##### comparison between my conditioning 1 and DDD conditioning.
# the difference is not coming from the division with n, because if lambda2=mu2=0 that difference should not count, but it does.
brtsM      <- c(-10, -10)
brtsS      <- c(-4)
lambdas    <- c(0.5, 0)
mus        <- c(0.4, 0)

Pc <- sls:::Pc_1shift(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus); Pc$Pc1
Pc1_DDD <- sls:::DDD_conditioning(brtsM = brtsM, brtsS = brtsS, lambdas = lambdas, mus = mus, lx = 300); Pc1_DDD #ddep=1 and ddep=3 are the same for constant rates

#####
#my function
PM1
PS
PM2cs
PM2cp
