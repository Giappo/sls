#' @title BISSE's D
#' @author Giovanni Laudanno
#' @description Provides BISSE's D function
#' @inheritParams default_params_doc
#' @param pars parameters of the likelihood functions:
#' \itemize{
#'   \item pars[1] is lambda, i.e. speciation rate;
#'   \item pars[2] is mu, i.e. extinction rate;
#' }
#' @return D(t)
#' @export
Dt <- function(pars, t_0, t_f, E0, D0) {
  lambda <- pars[1]
  mu     <- pars[2]
  TT <- t_f - t_0
  LL <- exp(
    (mu - lambda) * TT
  )
  FF <- lambda * (1 - E0) - LL * (mu - E0 * lambda)
  DD <- (LL * D0 * (lambda - mu) ^ 2) / FF ^ 2
  return(DD)
}

#' @title BISSE's E
#' @author Giovanni Laudanno
#' @description Provides BISSE E's function
#' @inheritParams default_params_doc
#' @param pars parameters of the likelihood functions:
#' \itemize{
#'   \item pars[1] is lambda, i.e. speciation rate;
#'   \item pars[2] is mu, i.e. extinction rate;
#' }
#' @return E(t)
#' @export
Et <- function(pars, t_0, t_f, E0, D0) {
  lambda <- pars[1]
  mu     <- pars[2]
  TT <- t_f - t_0
  LL <- exp(
    (mu - lambda) * TT
  )
  FF <- lambda * (1 - E0) - LL * (mu - E0 * lambda)
  EE <- 1 - (1 - E0) * (lambda - mu) / FF
  return(EE)
}

#' @title BISSE loglik
#' @author Giovanni Laudanno
#' @description Provides BISSE loglik function
#' @inheritParams default_params_doc
#' @param pars parameters of the likelihood functions:
#' \itemize{
#'   \item pars[1] is lambda, i.e. speciation rate;
#'   \item pars[2] is mu, i.e. extinction rate;
#' }
#' @return loglik
#' @export
loglik_bisse <- function(
  pars,
  brts,
  n_0 = 2,
  t_ds = NULL,
  D0s = NULL,
  t_p = 0
) {

  testit::assert(length(t_ds) == length(D0s))
  testit::assert(all(brts > t_p))

  lambda <- pars[1]

  tips <- (n_0 - 1) + length(brts) - length(t_ds)

  BRTS <- sort(c(brts, t_p, t_ds), decreasing = TRUE)
  maxt <- length(BRTS); mint <- 2; times <- maxt:mint; times
  lefts  <- rep(1, length(times))
  rights <- rep(2, length(times))
  for (t in times) {
    if (t == maxt) {
      D0 <- rep(1, tips)
      E0 <- 0
    } else {
      D0 <- DD
      E0 <- EE
    }
    l_d   <- length(D0)
    pool <- 1:l_d
    DD   <- rep(NA, l_d)
    for (N in pool) {
      t_0 <- BRTS[t]; t_f <- BRTS[t - 1]
      DD[N] <- Dt(pars = pars, t_0 = t_0, t_f = t_f, E0 = E0, D0 = D0[N])
    }
    EE    <- Et(pars = pars, t_0 = t_0, t_f = t_f, E0 = E0, D0 = D0)
    left  <- lefts[t - 1]; right <- rights[t - 1]

    if (t_f %in% t_ds) {
      if (length(t_ds) == 1) {
        DS0 <- D0s
      } else {
        DS0 <- D0s[which[t_ds == t_f]]
      }
      DD <- c(DD, DS0)
    } else {
      if (length(DD) > 1) {
        DD <- c(
          lambda ^ (t != mint) * DD[left] * DD[right],
          DD[-c(left, right)]
        )
      }
    }

  }
  DDf <- prod(DD)
  return(log(DDf))
}

#' @title BISSE loglik
#' @author Giovanni Laudanno
#' @description Provides BISSE loglik function (alternative version)
#' @inheritParams default_params_doc
#' @return loglik
#' @export
loglik_bisse2 <- function(
  pars,
  brts,
  n_0 = 2,
  t_0 = 0,
  E0 = 0,
  D0 = 1,
  LOG = TRUE,
  lambdaterms = TRUE
) {
  lambda <- pars[1]
  BRTS <- c(rep(brts[1], n_0 - 1), brts)
  DD <- prod(Dt(pars = pars, t_f = BRTS, t_0 = t_0, E0 = E0, D0 = D0))
  DD <- DD * lambda ^ (length(brts[-1]) * lambdaterms)
  out <- (LOG) * log(DD) + (1 - LOG) * DD
  return(out)
}
