#' @title Fourier term
#' @author Giovanni Laudanno
#' @description Provides fourier terms: exp((2 * pi * i)/N)
#' @inheritParams default_params_doc
#' @return Fourier term
#' @export
phase.factor <- function(N, n, k) {exp(1i * 2 * pi * n * k * (N^-1))}

#' @title DFT
#' @author Giovanni Laudanno
#' @description Computes the Discrete Fourier Transform (DFT)
#' @inheritParams default_params_doc
#' @return DFT
#' @export
DFT <- function(vec) {
  if (is.matrix(vec)) {Nmax <- nrow(vec)}
  if (is.vector(vec)) {Nmax <- length(vec)}
  O <- outer(X = 1:Nmax, Y = 1:Nmax, FUN = function(n,k) sls::phase.factor(n = n, k = k, N = Nmax))
  O %*% vec
}

#' @title Inverse DFT
#' @author Giovanni Laudanno
#' @description Computes the inverse Discrete Fourier Transform
#' @inheritParams default_params_doc
#' @return Inverse DFT
#' @export
IDFT <- function(vec) {
  if (is.matrix(vec)) {Nmax <- nrow(vec)}
  if (is.vector(vec)) {Nmax <- length(vec)}
  O <- outer(X = 1:Nmax, Y = 1:Nmax, FUN = function(n,k) sls::phase.factor(n = n, k = k, N = Nmax))
  solve(O) %*% vec
}

#' @title Combine pn
#' @author Giovanni Laudanno
#' @description Convolutes all the processes before the shift and imposes the death before the present of all species that are not visible in the phylogeny
#' @inheritParams default_params_doc
#' @return Convolution of the probabilities for all the processes
#' @export
combine_pns <- function(lambda, mu, ts, tbar, nmax = 1e2, fun = sls:::pn_bar){
  nvec <- 1:nmax
  N <- length(ts)
  X <- vector("list", N)
  for (t in 1:N)
  {
    X[[t]] <- fun(n = nvec, t = ts[t], lambda = lambda, mu = mu, tbar = tbar)
  }
  pippo <- matrix(unlist(lapply(X, FUN = sls::DFT)), nrow = N, byrow = T); rownames(pippo) <- paste0("t", 1:N)
  # apply(pippo, MARGIN = 2, "prod")
  Re(sum((nvec^-1) * sls::IDFT(apply(pippo, MARGIN = 2, "prod")))) #awesome!
  # sum((nvec^-1)*(Re(IDFT(DFT(x0) * DFT(y0))))) #cool!
}

#' @title Combine pn (basic)
#' @author Giovanni Laudanno
#' @description Convolution function. Basic but slow.
#' @inheritParams default_params_doc
#' @return Convolution of the probabilities for all the processes
#' @export
combine_pns0 <- function(lambda, mu, ts, tbar, nmax = 1e2) {
  N  <- length(ts)

  ls <- rep(lambda, N)
  ms <- rep(mu, N)
  nvec <- 1:nmax

  ns      <- expand.grid(replicate(expr = nvec, n = N, simplify = FALSE)); colnames(ns) <- paste0("n",1:N); head(ns)
  LAMBDAS <- matrix(ls, nrow = dim(ns)[1], ncol = dim(ns)[2], byrow = T); colnames(LAMBDAS) <- paste0("lambda",1:N); head(LAMBDAS)
  MUS     <- matrix(ms, nrow = dim(ns)[1], ncol = dim(ns)[2], byrow = T); colnames(MUS) <- paste0("mu",1:N); head(MUS)
  TS      <- matrix(ts, nrow = dim(ns)[1], ncol = dim(ns)[2], byrow = T); colnames(TS) <- paste0("t",1:N); head(TS)

  out <- sum(
    apply(sls:::pn_bar(n = ns, t = TS, tbar = tbar, lambda = LAMBDAS, mu = MUS), MARGIN = 1, FUN = prod) *
      apply(ns, MARGIN = 1, FUN = sum)^-1
  ); out
  return(out)
}

#' @title Combine pn in the old wrong way
#' @author Giovanni Laudanno
#' @description Convolutes all the processes before the shift and imposes the death before the present of all species that are not visible in the phylogeny. It doesn't divide by N. Used to check on old models.
#' @inheritParams default_params_doc
#' @return Convolution of the probabilities for all the processes
#' @export
combine_pns_nodivision <- function(lambda, mu, ts, tbar, nmax = 1e2, fun = sls:::pn_bar){
  nvec <- 1:nmax
  N <- length(ts)
  X <- vector("list", N)
  for (t in 1:N)
  {
    X[[t]] <- fun(n = nvec, t = ts[t], lambda = lambda, mu = mu, tbar = tbar)
  }
  pippo <- matrix(unlist(lapply(X, FUN = sls::DFT)), nrow = N, byrow = T); rownames(pippo) <- paste0("t", 1:N)
  # apply(pippo, MARGIN = 2, "prod")
  Re(sum(sls::IDFT(apply(pippo, MARGIN = 2, "prod")))) #awesome!
  # sum((nvec^-1)*(Re(IDFT(DFT(x0) * DFT(y0))))) #cool!
}



#' #' Does something
#' #' @inheritParams default_params_doc
#' #' @return Convolution of the probabilities for all the processes
#' #' @export
#' combine_pns2 <- function(lambda, mu, ts, tbar, nmax = 1e2, fun = sls:::pn_bar){
#'   nvec <- 1:nmax
#'   N <- length(ts)
#'   X <- vector("list", N)
#'   for (t in 1:N) {
#'     X[[t]] <- fun(n = nvec, t = ts[t], lambda = lambda, mu = mu, tbar = tbar)
#'   }
#'   pippo <- matrix(unlist(lapply(X, FUN = fft)), nrow = N, byrow = T); rownames(pippo) <- paste0("t",1:N)
#'   # apply(pippo, MARGIN = 2, "prod")
#'   Re(sum((nvec^-1) * fft(apply(pippo, MARGIN = 2, "prod"), inverse = TRUE))) #awesome!
#'   # sum((nvec^-1)*(Re(IDFT(DFT(x0) * DFT(y0))))) #cool!
#' }

