#' @title Pt
#' @author Giovanni Laudanno
#' @description Nee's function: pt
#' @inheritParams default_params_doc
#' @return pt
#' @export
pt  <- function(lambda, mu, t) {
  time <- t
  exp_term <- exp(
    (mu - lambda) * time
  )
  out    <- (lambda == mu) * (1 / (1 + lambda * time)) +
    (lambda != mu) * (
      (lambda - mu + (lambda == mu)) /
        (lambda - mu * exp_term * (lambda != mu) + (lambda == mu))
    )
  return(unname(out))
}

#' @title ut
#' @author Giovanni Laudanno
#' @description Nee's function: ut
#' @inheritParams default_params_doc
#' @return ut
#' @export
ut  <- function(lambda, mu, t) {
  time <- t
  exp_term <- exp(
    (mu - lambda) * time
  )
  out    <- (lambda == mu) * (lambda * time / (1 + lambda * time)) +
    (lambda != mu) * (
      (lambda - lambda * exp_term + (lambda == mu)) /
        (lambda - mu * exp_term * (lambda != mu) + (lambda == mu))
    )
  return(unname(out))
}

#' @title Pn
#' @author Giovanni Laudanno
#' @description Nee's function: pn
#' @inheritParams default_params_doc
#' @return pn
#' @export
pn <- function(lambda, mu, t, n) {
  out <- (n > 0) * sls::pt(t = t, lambda = lambda, mu = mu) *
    (1 - sls::ut(t = t, lambda = lambda, mu = mu)) *
    sls::ut(t = t, lambda = lambda, mu = mu) ^ (n - 1 + 2 * (n == 0)) +
    (n == 0) * (1 - sls::pt(t = t, lambda = lambda, mu = mu))
  return(out)
}

#' @title Pn accounting for extinctions after the shifts
#' @author Giovanni Laudanno
#' @description Combine pn from Nee et al. and imposes the extinction before the present of all species not visible in the phylogeny
#' @inheritParams default_params_doc
#' @return pn times probability of extinction for n-1 species after the shift
#' @export
pn_bar <- function(lambda, mu, t, n, tbar = 0) {
  out <- (n > 0) * sls::pt(t = t, lambda = lambda, mu = mu) *
    (1 - sls::ut(t = t, lambda = lambda, mu = mu)) *
    n *
    sls::ut(t = t, lambda = lambda, mu = mu) ^ (n - 1) *
    (1 - sls::pt(t = tbar, lambda = lambda, mu = mu)) ^ (n - 1 + (n == 0)) +
    (n == 0) * (1 - sls::pt(t = t, lambda = lambda, mu = mu))
  return(out)
}
