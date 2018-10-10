Dt <- function(pars, t0, tf, E0, D0) {
  lambda <- pars[1]
  mu     <- pars[2]
  TT <- tf - t0
  LL <- exp((mu - lambda) * TT)
  FF <- lambda * (1 - E0) - LL * (mu - E0 * lambda)
  DD <- (LL * D0 * (lambda - mu)^2)/FF^2
  return(DD)
}

Et <- function(pars, t0, tf, E0, D0) {
  lambda <- pars[1]
  mu     <- pars[2]
  TT <- tf - t0
  LL <- exp((mu - lambda) * TT)
  FF <- lambda * (1 - E0) - LL * (mu - E0 * lambda)
  EE <- 1 - ((1 - E0) * (lambda - mu))/FF
  return(EE)
}

BISSE_loglik <- function(pars, brts, N0 = 2) {

  lambda <- pars[1]
  mu     <- pars[2]

  kvec <- (N0 - 1) + cumsum(brts == brts)
  tips <- kvec[length(kvec)]

  BRTS <- c(brts, 0)
  maxt <- length(BRTS); mint <- 2; times <- maxt:mint
  lefts  <- rep(1, length(times))
  rights <- rep(2, length(times))
  for (t in maxt:mint)
  {
    if (t == maxt)
    {
      D0 <- rep(1, tips)
      E0 <- 0
    }else
    {
      D0 <- DD
      E0 <- EE
    }
    lD   <- length(D0)
    pool <- 1:lD
    DD   <- rep(NA, lD)
    for (N in pool)
    {
     t0 <- BRTS[t]; tf <- BRTS[t - 1]
     DD[N] <- Dt(pars = pars, t0 = t0, tf = tf, E0 = E0, D0 = D0[N])
    }
    EE    <- Et(pars = pars, t0 = t0, tf = tf, E0 = E0, D0 = D0)
    left  <- lefts[t - 1]; right <- rights[t - 1]
    DD    <- c(lambda^(t != mint) * DD[left] * DD[right], DD[-c(left, right)])
  }; DD
  return(log(DD))
}

BISSE_loglik_M <- function(pars, brts, t_d, DS0 = 1) {

  lambda <- pars[1]
  mu     <- pars[2]

  # kvec <- 1 + cumsum(brts == brts)
  # tips <- kvec[length(kvec)] - 1
  tips <- 1 + length(brts) - length(t_d)

  BRTS <- sort(c(brts, 0, t_d), decreasing = TRUE)
  maxt <- length(BRTS); mint <- 2; times <- maxt:mint
  lefts  <- rep(1, length(times))
  rights <- rep(2, length(times))
  for (t in maxt:mint)
  {
    if (t == maxt)
    {
      D0 <- rep(1, tips)
      E0 <- 0
    }else
    {
      D0 <- DD
      E0 <- EE
    }
    lD   <- length(D0)
    pool <- 1:lD
    DD   <- rep(NA, lD)
    t0 <- BRTS[t]; tf <- BRTS[t - 1]
    for (N in pool)
    {
      DD[N] <- Dt(pars = pars, t0 = t0, tf = tf, E0 = E0, D0 = D0[N])
    }
    EE    <- Et(pars = pars, t0 = t0, tf = tf, E0 = E0, D0 = D0)
    left  <- lefts[t - 1]; right <- rights[t - 1]

    if (tf == t_d)
    {
      DD <- c(DD, DS0)
    }else
    {
      DD <- c(lambda^(t != mint) * DD[left] * DD[right], DD[-c(left, right)])
    }

  }; DD
  return(log(DD))
}

BISSE_loglik_S <- function(pars, brts) {

  lambda <- pars[1]
  mu     <- pars[2]

  if (length(brts) > 1)
  {
    log_D2 <- BISSE_loglik(pars = pars, brts = brts[-1]) + log(lambda)
    E2 <- Et(pars = pars, t0 = 0, tf = brts[2], E0 = 0, D0 = 1)
    DD <- Dt(pars = pars, t0 = brts[2], tf = brts[1], E0 = E2, D0 = exp(log_D2))
  }else{
    DD <- Dt(pars = pars, t0 = 0, tf = brts[1], E0 = 0, D0 = 1)
  }

  return(log(DD))
}

BISSE_loglik_shift_old <- function(parsM,
                               parsS,
                               brtsM,
                               brtsS) {
  loglikS <- sls::BISSE_loglik(pars = parsS, brts = brtsS, N0 = 1)
  loglik <- BISSE_loglik_M(pars = parsM, brts = brtsM, t_d = brtsS[1], DS0 = exp(loglikS))

  return(loglik)
}

BISSE_loglik_shift <- function(parsM,
                               parsS,
                               brtsM,
                               brtsS) {
  loglikS <- sls::BISSE_loglik(pars = parsS, brts = brtsS, N0 = 1)
  loglik  <- sls::BISSE_loglik(pars = parsM, brts = brtsM, N0 = 2,
                               tds = brtsS[1], D0s = exp(loglikS))

  return(loglik)
}


brtsM <- c(10, 8, 7, 4, 2)
brtsS <- c(3, 1, 0.5)
parsM <- c(0.3, 0.1)
parsS <- c(0.6, 0.05)
test1 <- BISSE_loglik_shift(parsM = parsM, parsS = parsS, brtsM = brtsM, brtsS = brtsS); test1
test2 <- DDD::dd_KI_loglik(pars1 = c(parsM, Inf, parsS, Inf, brtsS[1]),
                           pars2 = c(200, 1, 0, min(abs(brtsM[abs(brtsM) > pars1[7]])), 0, 2),
                           brtsM = brtsM,
                           brtsS = brtsS[-1],
                           missnumspec = c(0,0)
                           ); test2

diff2 <- function(parsM, parsS, brtsM, brtsS) {

  parsM1 <- parsM  ; parsS1 <- parsS;
  parsM2 <- parsM/2; parsS2 <- parsS*3/4;

  BISSE1 <- BISSE_loglik_shift(parsM = parsM1, parsS = parsS1, brtsM = brtsM, brtsS = brtsS); BISSE1
  BISSE2 <- BISSE_loglik_shift(parsM = parsM2, parsS = parsS2, brtsM = brtsM, brtsS = brtsS); BISSE2
  DDD1   <- DDD::dd_KI_loglik(pars1 = c(parsM1, Inf, parsS1, Inf, brtsS[1]),
                              pars2 = c(200, 1, 0, min(abs(brtsM[abs(brtsM) > pars1[7]])), 0, 2),
                              brtsM = brtsM,
                              brtsS = brtsS[-1],
                              missnumspec = c(0,0)
  ); DDD1
  DDD2   <- DDD::dd_KI_loglik(pars1 = c(parsM2, Inf, parsS2, Inf, brtsS[1]),
                              pars2 = c(200, 1, 0, min(abs(brtsM[abs(brtsM) > pars1[7]])), 0, 2),
                              brtsM = brtsM,
                              brtsS = brtsS[-1],
                              missnumspec = c(0,0)
  ); DDD2

  DeltaDDD   <- DDD1 - DDD2
  DeltaBISSE <- BISSE1 - BISSE2

  diff <- abs(DeltaBISSE - DeltaDDD)

  return(diff)
}

diff2(parsM = parsM, parsS = parsS, brtsM = brtsM, brtsS = brtsS)

BISSE_loglik_S(pars = pars, brts = brts)
sls::BISSE_loglik(pars = pars, brts = brts, N0 = 1)

BISSE_loglik_M(pars = parsM, brts = brtsM, t_d = brtsS[1], DS0 = 0.5)
sls::BISSE_loglik(pars = parsM, brts = brtsM, N0 = 2, tds = brtsS[1], D0s = 0.5)
