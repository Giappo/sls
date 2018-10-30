#' @title P-likelihood (with no division)
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convoluting Nee's functions. There is no division. It should yield the same likelihood as DDD.
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_slsP_nodiv <- function(parsM,
                              parsS,
                              brtsM,
                              brtsS,
                              cond,
                              N0 = 2,
                              nmax = 1e2
)
{

  lambdas <- c(parsM[1], parsS[1])
  mus     <- c(parsM[2], parsS[2])

  brtsM1 <- sort(abs(brtsM), decreasing = TRUE)
  brtsS1 <- sort(abs(brtsS), decreasing = TRUE)
  td <- brtsS1[1]

  testit::assert(all(sign(brtsM1) == sign(brtsS1[1])))
  testit::assert(
    all(sign(brtsS1 * !is.null(brtsS1)) == sign(td * !is.null(brtsS1)))
  )

  BRTSM <- rbind(brtsM1, rep(1, length(brtsM1))); dim(BRTSM) <- c(2, length(brtsM1))
  TD <- c(td, -1); dim(TD) <- c(2, 1)
  EVENTSM <- (M <- cbind(BRTSM, TD))[,order(-M[1,])]
  kvecM_after <- (N0 - 1) + cumsum(EVENTSM[2,])
  kvecM_before <- c(N0 - 1, kvecM_after[-length(kvecM_after)])
  kvecM_before1 <- kvecM_before[EVENTSM[2,] > 0]

  if (N0 == 2)
  {
    brtsM2 <- c(brtsM1[1], brtsM1)
  }else
  {
    brtsM2 <- brtsM1
  }
  tsM_pre_shift  <- brtsM2[brtsM2 > td] - td; tsM_pre_shift
  tsM_post_shift <- brtsM2[brtsM2 < td]     ; tsM_post_shift
  if (length(tsM_post_shift) == 0) {tsM_post_shift <- 0}
  if (length(tsM_pre_shift ) == 0) {cat("There are no branching times before the shift"); return(-Inf)}

  likM_pre_shift  <- sls::combine_pns_nodiv(
    lambda = lambdas[1],
    mu = mus[1],
    ts = tsM_pre_shift,
    tbar = td,
    nmax = nmax
  ); log(likM_pre_shift)
  likM_post_shift <- prod(
    sls::pn(n = 1,
            lambda = lambdas[1],
            mu = mus[1],
            t = tsM_post_shift
    )
  ) * sls:::pn(
    n = 1,
    t = td,
    lambda = lambdas[1],
    mu = mus[1]
  )^(length(tsM_pre_shift) - 1); log(likM_post_shift)
  likS_post_shift <- prod(
    sls::pn(
      n = 1,
      lambda = lambdas[2],
      mu = mus[2],
      t = brtsS1
    )
  ); log(likS_post_shift)
  loglikM0 <- log(likM_pre_shift) + log(likM_post_shift)
  loglikS0 <- log(likS_post_shift)


  # logcombinatoricsM <- ifelse(length(brtsM) > 1, sum(log(kvecM_before1[-1])), 0); exp(logcombinatoricsM)
  # logcombinatoricsS <- ifelse(length(brtsS) > 0, lfactorial(length(brtsS[-1]))   , 0); exp(logcombinatoricsS)
  logcombinatoricsM <- logcombinatoricsS <- 0 #combinatorics
  lM <- length(brtsM1[brtsM1 != brtsM1[1]]) # number of speciations in the Main clade
  lS <- length(brtsS1[brtsS1 != brtsS1[1]]) # number of speciations in the Sub clade
  loglikM <- loglikM0 + logcombinatoricsM + log(lambdas[1] + (length(brtsM) == 0)) * lM
  loglikS <- loglikS0 + logcombinatoricsS + log(lambdas[2] + (length(brtsS) == 0)) * lS

  Pc <- sls::Pc_1shift(
    parsM = parsM,
    parsS = parsS,
    brtsM = brtsM,
    brtsS = brtsS,
    cond = cond,
    nmax = nmax,
    N0 = N0
  )

  loglik <- loglikM + loglikS - log(Pc); loglik
  return(loglik)
}

#' @title Q-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood integrating the Q-equation. There is no division.
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_slsQ_nodiv <- function(parsM,
                        parsS,
                        brtsM,
                        brtsS,
                        cond,
                        N0 = 2,
                        nmax = 1e2
)
{

  lambdas <- c(parsM[1], parsS[1])
  mus     <- c(parsM[2], parsS[2])
  Ks      <- c(Inf, Inf)

  brtsM1 <- sort(abs(brtsM), decreasing = TRUE)
  brtsS1 <- sort(abs(brtsS), decreasing = TRUE)
  td <- brtsS1[1]

  missnumspec <- c(0,0)
  N0s <- c(N0, 1)

  testit::assert(all(sign(brtsM) == sign(td)))
  testit::assert(
    all(sign(brtsS * !is.null(brtsS)) == sign(td * !is.null(brtsS)))
  )

  #BASIC SETTINGS AND CHECKS
  Nclades <- length(lambdas)
  brtsM1 <- sort(c(0, abs(c(brtsM, td))), decreasing = TRUE)
  brtsS1 <- sort(c(0, abs(c(brtsS))), decreasing = TRUE)
  brts_list <- list(brtsM = brtsM1, brtsS = brtsS1)
  abstol <- 1e-16; reltol <- 1e-10
  nvec <- 0:nmax
  logliks <- rep(NA, Nclades)

  #LIKELIHOOD INTEGRATION
  clade <- 0 #clade == 1 is the main clade, clade == 2 is the subclade
  while ((clade <- clade + 1) <= Nclades)
  {
    #SETTING CLADE CONDITIONS
    lambda <- lambdas[clade]
    mu     <- mus[clade]
    K      <- Ks[clade]
    soc     <- N0s[clade]
    maxT   <- length(brts_list[[clade]])
    brts   <- brts_list[[clade]]

    #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
    Qi <- c(1, rep(0, nmax))
    Qt <- matrix(0, ncol = (nmax + 1), nrow = maxT)
    Qt[1,] <- Qi
    dimnames(Qt)[[2]] <- paste0("Q", 0:nmax)
    k <- soc
    t <- 2
    D <- C <- rep(1, maxT)

    #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
    while (t <= maxT)
    {
      #Applying A operator
      if (lambda == 0 && mu == 0)
      {
        Qt[t,] <- Qt[(t-1),]
      }else
      {
        transition_matrix <- DDD:::dd_loglik_M_aux(pars = c(lambda, mu, K), lx = nmax + 1, k = k, ddep = 1)
        Qt[t,] <- abs(expoRkit::expv(v = Qt[(t-1),], x = transition_matrix, t = abs(brts[t] - brts[t - 1])))
      }

      #Applying C operator (this is a trick to avoid precision issues)
      C[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * C[t]

      #Applying B operator
      if (t < maxT)
      {
        if (brts[t] != td)
        {
          # print("+1")
          Qt[t,] <- Qt[t,] * lambda
          # Qt[t,] <- Qt[t,] * k * lambda #combinatorics
          k <- k + 1
        }else
        {
          # print("-1")
          Qt[t,] <- Qt[t,] #* (k + nvec)^-1 #NODIV
          k <- k - 1
        }

        #Applying D operator (this works exactly like C)
        D[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * D[t]

        #Updating running parameters
        t <- t + 1
      }else{break}
    }

    #Selecting the state I am interested in
    vm <- 1/choose((k + missnumspec[clade]), k)
    P  <- vm * Qt[t, (missnumspec[clade] + 1)] #I have to include +1 because of Q0

    #Removing C and D effects from the LL
    loglik <- log(P) - sum(log(C)) - sum(log(D))

    #Various checks
    loglik <- as.numeric(loglik)
    if (is.nan(loglik) | is.na(loglik))
    {
      loglik <- -Inf
    }
    logliks[clade] <- loglik
  }

  Pc <- sls::Pc_1shift(
    parsM = parsM,
    parsS = parsS,
    brtsM = brtsM,
    brtsS = brtsS,
    cond = cond,
    nmax = nmax,
    N0 = N0
  )

  total_loglik <- sum(logliks) - log(Pc)
  return(total_loglik)
}

#' @title DDD-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood calling the routine from the DDD package
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_DDD <- function(parsM,
                       parsS,
                       brtsM,
                       brtsS,
                       cond,
                       N0 = 2,
                       nmax = 1e2
)
{
  pars1 = c(parsM[1], parsM[2], Inf, parsS[1], parsS[2], Inf, brtsS[1])
  pars2 = c(nmax,  #maximum number of species involved in the computation
            1,  #ddmodel: not actually used by this function
            0,  #conditioning
            min(abs(brtsM[abs(brtsM) > pars1[7]])), # tshift
            0,  #print things: not actually used by this function
            N0) #stem or crown (soc)
  brtsM1 = brtsM
  brtsS1 = brtsS[-1]
  missnumspec = c(0,0)

  testit::assert(pars2[4] %in% abs(brtsM1)) #tsplit in brtsM

  loglik <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = abs(brtsM1),
    brtsS = abs(brtsS1),
    missnumspec = missnumspec
  ); loglik

  Pc <- sls::Pc_1shift(
    parsM = parsM,
    parsS = parsS,
    brtsM = brtsM,
    brtsS = brtsS,
    cond = cond,
    nmax = nmax,
    N0 = N0
  ); Pc

  loglikC <- loglik - log(Pc)
  return(loglikC)
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
                               cond,
                               N0 = 2,
                               nmax = 1e2) {

  loglikS <- sls::loglik_bisse(pars = parsS, brts = brtsS, N0 = 1)
  loglik  <- sls::loglik_bisse(pars = parsM, brts = brtsM, N0 = N0,
                               tds = brtsS[1], D0s = exp(loglikS))

  Pc <- sls::Pc_1shift(
    parsM = parsM,
    parsS = parsS,
    brtsM = brtsM,
    brtsS = brtsS,
    cond = cond,
    nmax = nmax,
    N0 = N0
  ); Pc

  loglikC <- loglik - log(Pc)
  return(loglikC)
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
