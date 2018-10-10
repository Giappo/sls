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
BISSE_loglik <- function(pars, brts, N0 = 2, tds = NULL, D0s = NULL) {

  testit::assert(length(tds) == length(D0s))

  lambda <- pars[1]
  mu     <- pars[2]

  # kvec <- (N0 - 1) + cumsum(brts == brts)
  # tips <- kvec[length(kvec)]

  tips <- (N0 - 1) + length(brts) - length(tds)

  BRTS <- sort(c(brts, 0, tds), decreasing = TRUE)
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
  return(log(DD))
}

#' @title BISSE loglik with shift
#' @author Giovanni Laudanno
#' @description Provides BISSE loglik function in a presence of a shift.
#' @inheritParams default_params_doc
#' @return loglik
#' @export
BISSE_loglik_shift <- function(parsM,
                               parsS,
                               brtsM,
                               brtsS) {
  loglikS <- sls::BISSE_loglik(pars = parsS, brts = brtsS, N0 = 1)
  loglik  <- sls::BISSE_loglik(pars = parsM, brts = brtsM, N0 = 2,
                               tds = brtsS[1], D0s = exp(loglikS))

  return(loglik)
}
