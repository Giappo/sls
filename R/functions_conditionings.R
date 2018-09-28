#' @title DDD-conditioning
#' @author Giovanni Laudanno
#' @description Calculates the conditioning from the DDD package for the Key Innovation model
#' @inheritParams default_params_doc
#' @return Conditioning probability
#' @export
DDD_conditioning <- function(pars, brtsM, brtsS, lx, ddep = 1) {

  lambdas <- c(pars[1], pars[3])
  mus     <- c(pars[2], pars[4])
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
  probs[1, 1:lx] = 0
  probs[1:lx, 1] = 0
  auxM1 = rep(0:(lx - 1), times = lx) + rep(0:(lx - 1), each = lx)
  probs = probs * rep(0:(lx - 1), times = lx)/auxM1
  dim(probs) = c(lx, lx)
  probs = rbind(probs[2:lx, 1:lx], rep(0, lx))
  dim(probs) = c(lx * lx, 1)
  probs = DDD:::dd_loglik_M2(pars = pars1[1:3], lx = lx,
                             ddep = ddep, tt = abs(tpres - tinn), p = probs)
  dim(probs) = c(lx, lx)
  PM12 = sum(probs[2:lx, 2:lx])
  PM2 = sum(probs[1, 2:lx])
  logliknorm = log(2) + log(PM12 + PS * PM2); exp(logliknorm)
  return(exp(logliknorm))
}

#' @title sls-conditioning (old version)
#' @author Giovanni Laudanno
#' @description Calculates three different kind of conditioning probabilities
#' @inheritParams default_params_doc
#' @return Three different conditioning probabilities
#' @export
Pc_1shift0 <- function(pars, brtsM, brtsS, Nmax = 100) {

  lambdas <- c(pars[1], pars[3])
  mus     <- c(pars[2], pars[4])

  if (length(brtsS) > 0)
  {
    tc <- brtsM[1]
    ts <- brtsS[1]
    tp <- 0
    A <- abs(ts - tc); B <- abs(tp - ts)
    PS <- 1 - pn(n = 0, t = B, lambda = lambdas[2], mu = mus[2])
    nvec <- 1:Nmax
    ns1  <- row(matrix(NA, nrow = Nmax, ncol = Nmax))
    ns2  <- col(matrix(NA, nrow = Nmax, ncol = Nmax))
    pA   <- sls::pt(t = A, lambda = lambdas[1], mu = mus[1]); pA
    uA   <- sls::ut(t = A, lambda = lambdas[1], mu = mus[1]); uA
    pB1  <- sls::pt(t = B, lambda = lambdas[1], mu = mus[1]); pB1
    pB2  <- sls::pt(t = B, lambda = lambdas[2], mu = mus[2]); pB2
    pns1 <- sls::pn(n = ns1, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns1) <- paste0("ns1=", nvec); colnames(pns1) <- paste0("ns2=", nvec); head(pns1)
    pns2 <- sls::pn(n = ns2, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns2) <- paste0("ns1=", nvec); colnames(pns2) <- paste0("ns2=", nvec); head(pns2)
    aux1 <- pns1 * pns2 * (ns1/(ns1 + ns2)) * (1 - (1 - pB1)^ns2)
    P1   <- sum(aux1) #branch 2 survives till the present
    aux2 <- aux1 * (1 - (1 - pB1)^(ns1 - 1)); head(aux2)
    P2   <- sum(aux2) #both branches 1 and 2 survive till the present
  }else
  {
    P1 <- (1 - sls::pn(n = 0, lambda = lambdas[1], mu = mus[1]))
    P2 <- (1/2) * (1 - sls::pn(n = 0, lambda = lambdas[1], mu = mus[1]))^2 #both branches 1 and 2 survive till the present
    PS <- 0
  }

  Pc1  <- 2 * PS * P1 + 2 * (1 - PS) * P2
  Pc2  <- 2 * PS * P2
  Pc3  <- 2 * PS * P1

  return(list(Pc1 = Pc1, Pc2 = Pc2, Pc3 = Pc3))
}

#' @title sls-conditioning
#' @author Giovanni Laudanno
#' @description Calculates three different kind of conditioning probabilities
#' @inheritParams default_params_doc
#' @return Three different conditioning probabilities
#' @export
Pc_1shift <- function(pars1, pars2, brtsM, brtsS) {

  missnumspec <- c(0,0)
  lambdas <- c(pars1[1], pars1[4])
  mus     <- c(pars1[2], pars1[5])
  Ks      <- c(pars1[3], pars1[6])
  td      <- pars1[7]
  td      <- abs(td) * sign(brtsM[1])

  lx     <- pars2[1]
  cond   <- pars2[3]
  tsplit <- pars2[4]
  soc    <- pars2[6]

  tc <- brtsM[1]
  ts <- td
  tp <- 0
  A <- abs(ts - tc); B <- abs(tp - ts)

  if (soc != 2) {stop("Pc can be calculated only if phylogeny starts with a crown!")}

  PS   <- 1 - pn(n = 0, t = B, lambda = lambdas[2], mu = mus[2])
  nvec <- 1:lx
  ns1  <- row(matrix(NA, nrow = lx, ncol = lx))
  ns2  <- col(matrix(NA, nrow = lx, ncol = lx))
  pA   <- sls::pt(t = A, lambda = lambdas[1], mu = mus[1]); pA
  uA   <- sls::ut(t = A, lambda = lambdas[1], mu = mus[1]); uA
  pB1  <- sls::pt(t = B, lambda = lambdas[1], mu = mus[1]); pB1
  pB2  <- sls::pt(t = B, lambda = lambdas[2], mu = mus[2]); pB2
  pns1 <- sls::pn(n = ns1, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns1) <- paste0("ns1=", nvec); colnames(pns1) <- paste0("ns2=", nvec); head(pns1)
  pns2 <- sls::pn(n = ns2, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns2) <- paste0("ns1=", nvec); colnames(pns2) <- paste0("ns2=", nvec); head(pns2)
  aux1 <- pns1 * pns2 * (ns1/(ns1 + ns2)) * (1 - (1 - pB1)^ns2)
  P1   <- sum(aux1) #branch 2 survives till the present
  aux2 <- aux1 * (1 - (1 - pB1)^(ns1 - 1)); head(aux2)
  P2   <- sum(aux2) #both branches 1 and 2 survive till the present

  Pc1  <- 2 * PS * P1 + 2 * (1 - PS) * P2
  Pc2  <- 2 * PS * P2
  Pc3  <- 2 * PS * P1

  return(list(Pc1 = Pc1, Pc2 = Pc2, Pc3 = Pc3))
}
