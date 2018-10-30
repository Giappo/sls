#' @title BISSE's D
#' @author Giovanni Laudanno
#' @description Provides BISSE's D function
#' @inheritParams default_params_doc
#' @return D(t)
#' @export
Dt <- function(pars, t0, tf, E0, D0) {
  lambda <- pars[1]
  mu     <- pars[2]
  TT <- tf - t0
  LL <- exp((mu - lambda) * TT)
  FF <- lambda * (1 - E0) - LL * (mu - E0 * lambda)
  DD <- (LL * D0 * (lambda - mu)^2)/FF^2
  return(DD)
}

#' @title BISSE's E
#' @author Giovanni Laudanno
#' @description Provides BISSE E's function
#' @inheritParams default_params_doc
#' @return E(t)
#' @export
Et <- function(pars, t0, tf, E0, D0) {
  lambda <- pars[1]
  mu     <- pars[2]
  TT <- tf - t0
  LL <- exp((mu - lambda) * TT)
  FF <- lambda * (1 - E0) - LL * (mu - E0 * lambda)
  EE <- 1 - ((1 - E0) * (lambda - mu))/FF
  return(EE)
}

#' @title BISSE loglik
#' @author Giovanni Laudanno
#' @description Provides BISSE loglik function
#' @inheritParams default_params_doc
#' @return loglik
#' @export
loglik_bisse <- function(pars, brts, N0 = 2, tds = NULL, D0s = NULL, tp = 0) {

  testit::assert(length(tds) == length(D0s))
  testit::assert(all(brts > tp))

  lambda <- pars[1]
  mu     <- pars[2]

  # kvec <- (N0 - 1) + cumsum(brts == brts)
  # tips <- kvec[length(kvec)]

  tips <- (N0 - 1) + length(brts) - length(tds)

  BRTS <- sort(c(brts, tp, tds), decreasing = TRUE)
  maxt <- length(BRTS); mint <- 2; times <- maxt:mint; times
  lefts  <- rep(1, length(times))
  rights <- rep(2, length(times))
  for (t in times)
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

    if (tf %in% tds)
    {
      if (length(tds) == 1)
      {
        DS0 <- D0s
      }else
      {
        DS0 <- D0s[which[tds == tf]]
      }
      DD <- c(DD, DS0)
    }else
    {
      if (length(DD) > 1)
      {
        DD <- c(lambda^(t != mint) * DD[left] * DD[right], DD[-c(left, right)]); #print(DD)
      }
    }


  }
  DDf <- prod(DD)
  return(log(DDf))
}

#' @title BISSE loglik with shift
#' @author Giovanni Laudanno
#' @description Provides BISSE loglik function in a presence of a shift.
#' @inheritParams default_params_doc
#' @return loglik
#' @export
loglik_bisse_shift <- function(parsM,
                               parsS,
                               brtsM,
                               brtsS,
                               N0M = 2) {
  loglikS <- sls::loglik_bisse(pars = parsS, brts = brtsS, N0 = 1)
  loglik  <- sls::loglik_bisse(pars = parsM, brts = brtsM, N0 = N0M,
                               tds = brtsS[1], D0s = exp(loglikS))

  return(loglik)
}

#' @title BISSE loglik
#' @author Giovanni Laudanno
#' @description Provides BISSE loglik function (alternative version)
#' @inheritParams default_params_doc
#' @return loglik
#' @export
loglik_bisse2 <- function(pars, brts, N0 = 2, t0 = 0,
                          E0 = 0, D0 = 1,
                          LOG = TRUE, lambdaterms = TRUE) {
  lambda <- pars[1]
  BRTS <- c(rep(brts[1], N0 - 1), brts)
  DD <- prod(Dt(pars = pars, tf = BRTS, t0 = t0, E0 = E0, D0 = D0))
  DD <- DD * lambda^(length(brts[-1]) * lambdaterms)
  out <- (LOG) * log(DD) + (1 - LOG) * DD
  return(out)
}

#' @title BISSE loglik shift
#' @author Giovanni Laudanno
#' @description Provides BISSE loglik shift function (alternative version). Yields the old (wrong) BISSE result for the Main Clade only.
#' @inheritParams default_params_doc
#' @return loglik
#' @export
loglik_bisse_shift2 <- function(pars, brts, N0 = 2, t0 = 0, td,
                                LOG = TRUE, lambdaterms = TRUE) {
  testit::assert(all(td != brts))
  lambda <- pars[1]
  brts1 <- brts[brts > td]; brts2 <- sort(c(td, brts[brts < td]), decreasing = TRUE)
  DD1 <- sls::loglik_bisse2(pars, brts1, N0 = N0, t0 = td,
                       E0 = sls::Et(pars = pars, t0 = t0, tf = td, E0 = 0, D0 = 1),
                       LOG = FALSE, lambdaterms = FALSE)
  DD2 <- sls::loglik_bisse2(pars, brts2, N0 = (N0 + length(brts1) - 1) - 1, t0 = t0,
                       LOG = FALSE, lambdaterms = FALSE)

  DD <- DD1 * DD2 * lambda^(length(brts[-1]) * lambdaterms)
  out <- (LOG) * log(DD) + (1 - LOG) * DD
  return(out)
}
