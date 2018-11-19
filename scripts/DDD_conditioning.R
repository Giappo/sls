DDD_conditioning <- function(lambdas, mus, brtsM, brtsS, lx) {
  pars1 <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf)
  ddep <- 1
  tcrown <- brtsM[1]
  tpres  <- 0
  tinn   <- brtsS[1]
  # lx = lxS
  lxS <- lxM <- lx
  probs = rep(0, lx)
  probs[2] = 1
  probs = DDD::dd_loglik_M(pars1[4:6], lx, 0, ddep,
                           tt = abs(tpres - tinn), probs)
  PS = 1 - probs[1]
  lx = lxM
  probs = matrix(0, lx, lx)
  probs[2, 2] = 1
  dim(probs) = c(lx * lx, 1)
  ly = lx^2
  probs = DDD::dd_loglik_M2(pars = pars1[1:3], lx = lx,
                            ddep = ddep, tt = abs(tinn - tcrown), p = probs)
  dim(probs) = c(lx, lx)
  probs[1, 1:lx] = 0
  probs[1:lx, 1] = 0
  auxM1 = rep(0:(lx - 1), times = lx) + rep(0:(lx - 1), each = lx)
  probs = probs * rep(0:(lx - 1), times = lx)/auxM1
  dim(probs) = c(lx, lx)
  probs = rbind(probs[2:lx, 1:lx], rep(0, lx))
  dim(probs) = c(lx * lx, 1)
  probs = DDD::dd_loglik_M2(pars = pars1[1:3], lx = lx,
                            ddep = ddep, tt = abs(tpres - tinn), p = probs)
  dim(probs) = c(lx, lx)
  PM12 = sum(probs[2:lx, 2:lx])
  PM2 = sum(probs[1, 2:lx])
  logliknorm = log(2) + log(PM12 + PS * PM2); exp(logliknorm)
  return(exp(logliknorm))
}
