#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
debug_DDD_conditioning <- function(lambdas, mus, brtsM, brtsS, lx, ddep = 1) {
  pars1 <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, brtsS[1])
  laM = pars1[1]
  muM = pars1[2]
  KM = pars1[3]
  laS = pars1[4]
  muS = pars1[5]
  KS = pars1[6]
  tinn = -abs(pars1[7])
  lmax = lx #pars2[1]
  # ddep = pars2[2]

  m = 0
  lxM = min(max(1 + m[1], 1 + ceiling(KM)), ceiling(lmax))
  lxS = min(max(1 + m[1], 1 + ceiling(KS)), ceiling(lmax))

  tcrown = brtsM[1]
  tpres  = 0
  # tinn <- brtsS[1]
  lx = lxS
  # lxS <- lxM <- lx
  probs = rep(0, lx)
  probs[2] = 1
  probs = DDD:::dd_loglik_M(pars1[4:6], lx, 0, ddep,
                            tt = abs(tpres - tinn), probs)
  PS = 1 - probs[1]
  lx = lxM
  probs = matrix(0, lx, lx)
  probs[2, 2] = 1
  dim(probs) = c(lx * lx, 1)
  ly = lx^2
  probs = DDD:::dd_loglik_M2(pars = pars1[1:3], lx = lx,
                             ddep = ddep, tt = abs(tinn - tcrown), p = probs)
  dim(probs) = c(lx, lx)

  PM2cs = sum(probs[1, 2:lx])

  probs[1, 1:lx] = 0
  probs[1:lx, 1] = 0
  auxM1 = rep(0:(lx - 1), times = lx) + rep(0:(lx - 1), each = lx)
  probs = probs * rep(0:(lx - 1), times = lx)/auxM1
  dim(probs) = c(lx, lx)
  probs = rbind(probs[2:lx, 1:lx], rep(0, lx))

  # PM2cs_after = sum(probs[1, 2:lx])

  dim(probs) = c(lx * lx, 1)
  probs = DDD:::dd_loglik_M2(pars = pars1[1:3], lx = lx,
                             ddep = ddep, tt = abs(tpres - tinn), p = probs)
  dim(probs) = c(lx, lx)
  PM12 = sum(probs[2:lx, 2:lx])
  PM2 = sum(probs[1, 2:lx])
  logliknorm = log(2) + log(PM12 + PS * PM2); exp(logliknorm)
  return(list(PS = PS, PM2 = PM2, PM12 = PM12, PM2cs = PM2cs))
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
debug_Pc_1shift <- function(brtsM, brtsS, lambdas, mus) {

  tp <- 0 ;tc <- brtsM[1]; ts <- brtsS[1]
  testit:::assert(tp > ts)
  testit:::assert(ts > tc)

  PM1   <- (1 - sls::pn(n = 0, t = tp - tc, lambda = lambdas[1], mu = mus[1]))
  PM2cs <- (1 - sls::pn(n = 0, t = ts - tc, lambda = lambdas[1], mu = mus[1]))
  PS    <- (1 - sls::pn(n = 0, t = tp - ts, lambda = lambdas[2], mu = mus[2]))

  P1cs <- sls::pn(n = 1, t = ts - tc, lambda = lambdas[1], mu = mus[1])
  P0sp <- sls::pn(n = 0, t = tp - ts, lambda = lambdas[1], mu = mus[1])
  ucs  <- sls::ut(t = ts - tc, lambda = lambdas[1], mu = mus[1])

  PM2cp <- P1cs *
    (1 - P0sp) *
    ucs *
    (1 - ucs)^-1 *
    (1 - ucs * P0sp)^-1

  Pc1 <- PM1 * (PS * PM2cs + (1 - PS) * PM2cp)
  Pc2 <- PM1 * PS * PM2cp
  Pc3 <- PM1 * PS * PM2cs
  return(list(PS = PS, PM1 = PM1, PM2cs = PM2cs, PM2cp = PM2cp))
}
