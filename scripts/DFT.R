install.packages("spatialfil")
library(spatialfil)
M <- array(1:100, c(10, 10)); M
Mfil <- applyFilter(x = M, kernel = convKernel(sigma = 1.4, k = 'gaussian'))

omega <- function(N,n,k) exp(1i * 2 * pi * n * k * (N^-1))
DFT <- function(vec) {
  if (is.matrix(vec)) {Nmax <- nrow(vec)}
  if (is.vector(vec)) {Nmax <- length(vec)}
  O <- outer(X = 1:Nmax, Y = 1:Nmax, FUN = function(n,k) omega(n = n, k = k, N = Nmax))
  O %*% vec
}
DFT_inv <- function(vec) {
  if (is.matrix(vec)) {Nmax <- nrow(vec)}
  if (is.vector(vec)) {Nmax <- length(vec)}
  O <- outer(X = 1:Nmax, Y = 1:Nmax, FUN = function(n,k) omega(n = n, k = k, N = Nmax))
  solve(O) %*% vec
}
x0 <- fun(n = nvec, t = ts[1], lambda = lambda, mu = mu, tbar = tbar)
y0 <- fun(n = nvec, t = ts[2], lambda = lambda, mu = mu, tbar = tbar)
sum(Re(DFT_inv(DFT(x0) * DFT(y0)))) #cool!

# n1n2 <- function(Nmax) {
#   n1 <- matrix(1:Nmax, nrow = Nmax, ncol = Nmax); rownames(n1) <- paste0("n1=",1:Nmax); colnames(n1) <- paste0("n2=",1:Nmax); n1
#   n2 <- t(n1); names(n2) <- names(n1); n2
#   Mn1n2 <- (n1 + n2)^-1; names(Mn1n2) <- names(n1); Mn1n2
# }
# nvec <- 1:30; Nmax <- length(nvec)
# x0 <- fun(n = nvec, t = ts[1], lambda = lambda, mu = mu, tbar = tbar)
# y0 <- fun(n = nvec, t = ts[2], lambda = lambda, mu = mu, tbar = tbar)
# x1 <- matrix(x0, nrow = Nmax, ncol = Nmax)
# y1 <- t(matrix(y0, nrow = Nmax, ncol = Nmax))
# sum(Re(DFT_inv(
#   DFT(x1 * n1n2(Nmax)) *
#   DFT(y1)
# ))) #awesome!
combine_pns <- function(lambda, mu, ts, tbar, nmax = 1e2, fun = sls:::pn_bar){
  nvec <- 1:nmax
  N <- length(ts)
  X <- vector("list", N)
  for (t in 1:N) {
    X[[t]] <- fun(n = nvec, t = ts[t], lambda = lambda, mu = mu, tbar = tbar)
  }
  pippo <- matrix(unlist(lapply(X, FUN = DFT)), nrow = N, byrow = T); rownames(pippo) <- paste0("t",1:N)
  # apply(pippo, MARGIN = 2, "prod")
  Re(sum((nvec^-1) * DFT_inv(apply(pippo, MARGIN = 2, "prod")))) #awesome!
  # sum((nvec^-1)*(Re(DFT_inv(DFT(x0) * DFT(y0))))) #cool!
}

lambda <- 0.3
mu <- 0.1
ts <- c(10, 8, 7, 6, 3)
tbar <- 2
nmax <- 15
combine_pns0(lambda = lambda, mu = mu, ts = ts, tbar = tbar, nmax = nmax)
combine_pns(lambda = lambda, mu = mu, ts = ts, tbar = tbar, nmax = nmax)
combine_pns2(lambda = lambda, mu = mu, ts = ts, tbar = tbar, nmax = nmax)

pn(lambda = lambda, mu = mu, t = ts[1], n = 1:5)
pn_bar(lambda = lambda, mu = mu, t = ts[1], n = 1:5)
