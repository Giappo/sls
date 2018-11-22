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
d_t <- function(pars, t_0, t_f, e_0, d_0) {
  lambda <- pars[1]
  mu     <- pars[2]
  delta_t <- t_f - t_0
  exp_factor <- exp(
    (mu - lambda) * delta_t
  )
  denominator <- lambda * (1 - e_0) - exp_factor * (mu - e_0 * lambda)
  DD <- (exp_factor * d_0 * (lambda - mu) ^ 2) / denominator ^ 2
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
e_t <- function(pars, t_0, t_f, e_0, d_0) {
  lambda <- pars[1]
  mu     <- pars[2]
  delta_t <- t_f - t_0
  exp_factor <- exp(
    (mu - lambda) * delta_t
  )
  denominator <- lambda * (1 - e_0) - exp_factor * (mu - e_0 * lambda)
  es_t <- 1 - (1 - e_0) * (lambda - mu) / denominator
  return(es_t)
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
  d_0s = NULL,
  t_p = 0
) {

  testit::assert(length(t_ds) == length(d_0s))
  testit::assert(all(brts > t_p))

  lambda <- pars[1]

  tips <- (n_0 - 1) + length(brts) - length(t_ds)

  BRTS <- sort(c(brts, t_p, t_ds), decreasing = TRUE)
  maxt <- length(BRTS); mint <- 2; times <- maxt:mint; times
  lefts  <- rep(1, length(times))
  rights <- rep(2, length(times))
  for (t in times) {
    if (t == maxt) {
      d_0 <- rep(1, tips)
      e_0 <- 0
    } else {
      d_0 <- ds_t
      e_0 <- es_t
    }
    l_d   <- length(d_0)
    pool <- 1:l_d
    ds_t   <- rep(NA, l_d)
    for (N in pool) {
      t_0 <- BRTS[t]; t_f <- BRTS[t - 1]
      ds_t[N] <- d_t(pars = pars, t_0 = t_0, t_f = t_f, e_0 = e_0, d_0 = d_0[N])
    }
    es_t    <- e_t(pars = pars, t_0 = t_0, t_f = t_f, e_0 = e_0, d_0 = d_0)
    left  <- lefts[t - 1]; right <- rights[t - 1]

    if (t_f %in% t_ds) {
      if (length(t_ds) == 1) {
        DS0 <- d_0s
      } else {
        DS0 <- d_0s[which[t_ds == t_f]]
      }
      ds_t <- c(ds_t, DS0)
    } else {
      if (length(ds_t) > 1) {
        ds_t <- c(
          lambda ^ (t != mint) * ds_t[left] * ds_t[right],
          ds_t[-c(left, right)]
        )
      }
    }

  }
  prod_d <- prod(ds_t)
  return(log(prod_d))
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
  e_0 = 0,
  d_0 = 1,
  log_scale = TRUE,
  lambdaterms = TRUE
) {
  lambda <- pars[1]
  BRTS <- c(rep(brts[1], n_0 - 1), brts)
  prod_d <- prod(d_t(pars = pars, t_f = BRTS, t_0 = t_0, e_0 = e_0, d_0 = d_0))
  prod_d <- prod_d * lambda ^ (length(brts[-1]) * lambdaterms)
  out <- (log_scale) * log(prod_d) + (1 - log_scale) * prod_d
  return(out)
}
