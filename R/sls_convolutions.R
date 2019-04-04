#' @title Fourier term
#' @author Giovanni Laudanno
#' @description Provides fourier terms: exp((2 * pi * i)/n_max)
#' @inheritParams default_params_doc
#' @return Fourier term
#' @export
phase_factor <- function(n_max, n, k) {
  exp(1i * 2 * pi * n * k * (n_max ^ -1))
}

#' @title dft
#' @author Giovanni Laudanno
#' @description Computes the Discrete Fourier Transform (dft)
#' @inheritParams default_params_doc
#' @return dft
#' @export
dft <- function(vec) {
  if (is.matrix(vec)) {
    n_max <- nrow(vec)
  }
  if (is.vector(vec)) {
    n_max <- length(vec)
  }
  oo <- outer(
    X = 1:n_max,
    Y = 1:n_max,
    FUN = function(n, k) sls::phase_factor(n = n, k = k, n_max = n_max)
  )
  oo %*% vec
}

#' @title Inverse dft
#' @author Giovanni Laudanno
#' @description Computes the inverse Discrete Fourier Transform
#' @inheritParams default_params_doc
#' @return Inverse dft
#' @export
idft <- function(vec) {
  if (is.matrix(vec)) {
    n_max <- nrow(vec)
  }
  if (is.vector(vec)) {
    n_max <- length(vec)
  }
  oo <- outer(
    X = 1:n_max,
    Y = 1:n_max,
    FUN = function(n, k) sls::phase_factor(n = n, k = k, n_max = n_max)
  )
  solve(oo) %*% vec
}

#' @title Combine pn
#' @author Giovanni Laudanno
#' @description Convolutes all the processes before the shift and
#'  imposes the death before the present of all species that
#'  are not visible in the phylogeny
#' @inheritParams default_params_doc
#' @return Convolution of the probabilities for all the processes
#' @export
combine_pns <- function(
  lambda,
  mu,
  times,
  tbar,
  n_max = 1e2,
  fun = sls::pn_bar
) {
  nvec <- 1:n_max
  n_t <- length(times)
  p_t_n <- vector("list", n_t)
  for (t in 1:n_t) {
    p_t_n[[t]] <- fun(
      n = nvec,
      t = times[t],
      lambda = lambda,
      mu = mu,
      tbar = tbar
    )
  }
  dft_p_t_n <- matrix(
    unlist(lapply(p_t_n, FUN = sls::dft)),
    nrow = n_t,
    byrow = TRUE
  )
  rownames(dft_p_t_n) <- paste0("t", 1:n_t)
  # I take the product of "dft_p_t_n" over all the branching times.
  # I transform back the result (which is a vector of length n_max)
  # using inverse dft. Finally I take the real part of the sum of
  # the back-transformed vector times the factor 1/n
  Re(sum(
    (nvec ^ -1) * sls::idft(apply(dft_p_t_n, MARGIN = 2, "prod"))
  )) # awesome!
}

#' @title Combine pn (basic)
#' @author Giovanni Laudanno
#' @description Convolution function. Basic but slow.
#' @inheritParams default_params_doc
#' @return Convolution of the probabilities for all the processes
#' @export
combine_pns0 <- function(
  lambda,
  mu,
  times,
  tbar,
  n_max = 1e2
) {
  n_t  <- length(times)

  ls <- rep(lambda, n_t)
  ms <- rep(mu, n_t)
  nvec <- 1:n_max

  ns      <- expand.grid(replicate(expr = nvec, n = n_t, simplify = FALSE))
  colnames(ns) <- paste0("n", 1:n_t)
  lambda_n <- matrix(ls, nrow = dim(ns)[1], ncol = dim(ns)[2], byrow = T)
  colnames(lambda_n) <- paste0("lambda", 1:n_t)
  mu_n     <- matrix(ms, nrow = dim(ns)[1], ncol = dim(ns)[2], byrow = T)
  colnames(mu_n) <- paste0("mu", 1:n_t)
  t_n   <- matrix(times, nrow = dim(ns)[1], ncol = dim(ns)[2], byrow = T)
  colnames(t_n) <- paste0("t", 1:n_t)

  out <- sum(
    apply(sls::pn_bar(
      n = ns,
      t = t_n,
      tbar = tbar,
      lambda = lambda_n,
      mu = mu_n
    ), MARGIN = 1, FUN = prod) *
      apply(ns, MARGIN = 1, FUN = sum) ^ -1
  ); out
  return(out)
}

#' @title Combine pn in the old wrong way
#' @author Giovanni Laudanno
#' @description Convolutes all the processes before the shift and
#'  imposes the death before the present of all species that
#'  are not visible in the phylogeny.
#'  It doesn't divide by N. Used to check on old models.
#' @inheritParams default_params_doc
#' @return Convolution of the probabilities for all the processes
#' @export
combine_pns_nodiv <- function(
  lambda,
  mu,
  times,
  tbar,
  n_max = 1e2,
  fun = sls::pn_bar
) {
  nvec <- 1:n_max
  n_t <- length(times)
  p_t_n <- vector("list", n_t)
  for (t in 1:n_t) {
    p_t_n[[t]] <- fun(
      n = nvec,
      t = times[t],
      lambda = lambda,
      mu = mu,
      tbar = tbar
    )
  }
  dft_p_t_n <- matrix(unlist(lapply(p_t_n, FUN = sls::dft)), nrow = n_t, byrow = T)
  rownames(dft_p_t_n) <- paste0("t", 1:n_t)
  Re(sum(sls::idft(apply(dft_p_t_n, MARGIN = 2, "prod")))) #awesome!
}
