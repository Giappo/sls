#analytical auxiliary function
ht <- function(pars, z = 0, t) {
  lambda <- pars[1]
  mu     <- pars[2]
  LL     <- exp((mu - lambda) * t)
  ht     <- (lambda * (1 - LL))/(lambda * (1 - z) - (mu - lambda * z) * LL)
  return(ht)
}

#alt for BISSE_loglik
BISSE_loglik2 <- function(pars, brts, N0, t0,
                          E0 = 0, D0 = 1,
                          LOG = TRUE, lambdaterms = TRUE) {
  lambda <- pars[1]
  BRTS <- c(rep(brts[1], N0 - 1), brts)
  DD <- prod(Dt(pars = pars, tf = BRTS, t0 = t0, E0 = E0, D0 = D0))
  DD <- DD * lambda^(length(brts[-1]) * lambdaterms)
  out <- (LOG) * log(DD) + (1 - LOG) * DD
  return(out)
}

#alt for BISSE_loglik with a split in the middle td
BISSE_loglik3 <- function(pars, brts, N0, t0, td,
                          LOG = TRUE, lambdaterms = TRUE) {
  testit::assert(all(td != brts))
  lambda <- pars[1]
  brts1 <- brts[brts > td]; brts2 <- sort(c(td, brts[brts < td]), decreasing = TRUE)
  DD1 <- BISSE_loglik2(pars, brts1, N0 = N0, t0 = td,
                       E0 = Et(pars = pars, t0 = t0, tf = td, E0 = 0, D0 = 1),
                       LOG = FALSE, lambdaterms = FALSE)
  DD2 <- BISSE_loglik2(pars, brts2, N0 = (N0 + length(brts1) - 1), t0 = 0,
                       LOG = FALSE, lambdaterms = FALSE)

  DD <- DD1 * DD2 * lambda^(length(brts[-1]) * lambdaterms)
  out <- (LOG) * log(DD) + (1 - LOG) * DD
  return(out)
}

#alt for BISSE_loglik with a decoupling in the middle td
BISSE_loglik4 <- function(pars, brts, N0, t0, td,
                          LOG = TRUE, lambdaterms = TRUE) {
  testit::assert(all(td != brts))
  lambda <- pars[1]
  brts1 <- brts[brts > td]; brts2 <- sort(c(td, brts[brts < td]), decreasing = TRUE)
  DD1 <- BISSE_loglik2(pars, brts1, N0 = N0, t0 = td,
                       E0 = Et(pars = pars, t0 = t0, tf = td, E0 = 0, D0 = 1),
                       LOG = FALSE, lambdaterms = FALSE)
  DD2 <- BISSE_loglik2(pars, brts2, N0 = (N0 + length(brts1) - 1) - 1, t0 = 0,
                       LOG = FALSE, lambdaterms = FALSE)

  DD <- DD1 * DD2 * lambda^(length(brts[-1]) * lambdaterms)
  out <- (LOG) * log(DD) + (1 - LOG) * DD
  return(out)
}

#BISSE version of slsP (?)
BISSE_loglik5 <- function(pars, brts, N0, t0, td,
                          LOG = TRUE, lambdaterms = TRUE) {
  testit::assert(all(td != brts))
  lambda <- pars[1]
  brts1 <- brts[brts > td]; brts2 <- sort(c(td, brts[brts < td]), decreasing = TRUE)
  DD1 <- BISSE_loglik2(pars, brts1, N0 = N0, t0 = td,
                       E0 = (Ed <- Et(pars = pars, t0 = t0, tf = td, E0 = 0, D0 = 1)),
                       LOG = FALSE, lambdaterms = FALSE)

  k <- (N0 + length(brts1) - 1)
  n <- 0:1000
  aux1.1 <- sum(
    choose(n + 2 * k + 1, n)  * (k + n)^-1 * (ht(pars = pars, t = td) * Ed)^n
  )
  aux1.2 <- sum(
    choose(n + 2 * k + 1, n)  * (ht(pars = pars, t = td) * Ed)^n
  )
  divterm1 <- aux1.1 / aux1.2

  aux2.1 <- sum(
    choose(n + 2 * k + 1, n)  * (k + n)^-1 * (ht(pars = pars, t = (brts[1] - td)) * Ed)^n
  )
  aux2.2 <- sum(
    choose(n + 2 * k + 1, n)  * (ht(pars = pars, t = (brts[1] - td)) * Ed)^n
  )
  divterm2 <- aux2.1 / aux2.2

  DD2 <- BISSE_loglik2(pars, brts2, N0 = k - 1, t0 = 0,
                       LOG = FALSE, lambdaterms = FALSE)

  DD <- DD1 * DD2 * lambda^(length(brts[-1]) * lambdaterms)
  out  <- (LOG) * log(DD) + (1 - LOG) * DD
  out1 <- out + log(divterm1)
  out2 <- out + log(divterm2)
  return(list(out1 = out1, out2 = out2))
}

##
pars = c(0.3, 0.2); brts = c(10, 6, 2); N0 = 2; t0 = 0
td <- 5; brts1 <- brts[brts > td]; brts2 <- c(td, brts[brts < td])
test1  <- BISSE_loglik(pars, brts); test1
test2  <- BISSE_loglik2(pars, brts, N0 = 2, t0 = t0, lambdaterms = T); test2

##
pars = c(0.3, 0.2); brts = c(10, 6, 2); N0 = 2; t0 = 0
test1 <- BISSE_loglik2(pars = pars, brts = brts, N0 = N0, t0 = t0); test1
test2 <- BISSE_loglik3(pars = pars, brts = brts, N0 = N0, t0 = t0, td = runif(n = 1, min = 1, max = 9)); test2

##
pars = c(0.3, 0.2); brts = c(10, 6, 2); N0 = 2; t0 = 0; td = runif(n = 1, min = 1, max = 9); td
test1 <- BISSE_loglik4(pars = pars, brts = brts, N0 = N0, t0 = t0, td = td); test1
test2 <- sls::BISSE_loglik_shift(parsM = pars, parsS = pars, brtsM = brts, brtsS = td, N0M = 2) -
  sls::BISSE_loglik(pars = pars, brts = td, N0 = 1); test2

##
pars = c(0.3, 0.2); brts = c(10, 6, 2); N0 = 2; t0 = 0; set.seed(1); td = runif(n = 1, min = 1, max = 9); td
test1 <- BISSE_loglik5(pars = pars, brts = brts, N0 = N0, t0 = t0, td = td); test1
test2 <- sls::loglik_slsP_nodivision(pars1 = c(parsM[1], parsM[2], Inf, parsM[1], parsM[2], Inf, td),
                                     pars2 = c(200, 1, 0, min(brts[brts > td]), 0, N0),
                                     missnumspec = 0, brtsM = brts, brtsS = NULL) -
         sls::BISSE_loglik(pars = pars, brts = td, N0 = 1); test2
