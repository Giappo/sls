#' @title P-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convoluting Nee's functions
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_slsP <- function(pars1, pars2 = c(100, 1, 1, brtsM[2], 0, 2), brtsM, brtsS, missnumspec = c(0,0)) {

  lambdas <- c(pars1[1], pars1[4])
  mus     <- c(pars1[2], pars1[5])
  Ks      <- c(pars1[3], pars1[6])
  td      <- pars1[7]

  nmax   <- pars2[1]
  cond   <- pars2[3]
  tsplit <- pars2[4]
  N0     <- pars2[6]

  dummy <- missnumspec #change this if you want to include missing species. it's probably not difficult, just modify pn for LM and LS

  brtsM1 <- sort(abs(brtsM), decreasing = TRUE)
  if (is.null(unlist(brtsS))) {brtsS1 <- td} else
  {
    brtsS1 <- sort(abs(c(brtsS, td)), decreasing = TRUE)
  }
  td <- abs(td) * sign(brtsM1[1])
  tsplit <- abs(tsplit) * sign(brtsM1[1])

  testit::assert(all(sign(brtsM1) == sign(td)))
  testit::assert(
    all(sign(brtsS1 * !is.null(brtsS1)) == sign(td * !is.null(brtsS1)))
  )
  testit::assert(all(sign(tsplit) == sign(td)))
  if (!(tsplit %in% brtsM1)) {stop('tsplit has to be in the main clade branching times!')}

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

  likM_pre_shift  <- sls::combine_pns(lambda = lambdas[1], mu = mus[1], ts = tsM_pre_shift, tbar = td, nmax = nmax); log(likM_pre_shift)
  likM_post_shift <- prod(
    sls::pn(n = 1, lambda = lambdas[1], mu = mus[1], t = tsM_post_shift)
  ) * sls:::pn(n = 1, t = td, lambda = lambdas[1], mu = mus[1])^(length(tsM_pre_shift) - 1); log(likM_post_shift)
  likS_post_shift <- prod(
    sls::pn(n = 1, lambda = lambdas[2], mu = mus[2], t = brtsS1        )
  ); log(likS_post_shift)
  loglikM0 <- log(likM_pre_shift) + log(likM_post_shift)
  loglikS0 <- log(likS_post_shift)

  logcombinatoricsM <- ifelse(length(brtsM) > 1, sum(log(kvecM_before1[-1])), 0)
  logcombinatoricsS <- ifelse(length(brtsS) > 0, lfactorial(length(brtsS))   , 0)
  lM <- length(brtsM1[brtsM1 != brtsM1[1]]) # number of speciations in the Main clade
  lS <- length(brtsS1[brtsS1 != brtsS1[1]]) # number of speciations in the Sub clade
  loglikM <- loglikM0 + logcombinatoricsM + log(lambdas[1] + (length(brtsM) == 0)) * lM
  loglikS <- loglikS0 + logcombinatoricsS + log(lambdas[2] + (length(brtsS) == 0)) * lS

  Pc <- 1
  if (!missing(cond))
  {
    if (cond %in% c(1,2,3))
    {
      conditioning <- sls::Pc_1shift(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS)
      Pc <- conditioning[[cond]]; Pc
    }
  }

  loglik <- loglikM + loglikS - log(Pc); loglik
  return(loglik)
}

#' @title P-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convoluting Nee's functions
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_slsP_nodivision <- function(pars1, pars2 = c(100, 1, 1, brtsM[2], 0, 2), brtsM, brtsS, missnumspec = c(0,0)) {

  lambdas <- c(pars1[1], pars1[4])
  mus     <- c(pars1[2], pars1[5])
  Ks      <- c(pars1[3], pars1[6])
  td      <- pars1[7]

  nmax   <- pars2[1]
  cond   <- pars2[3]
  tsplit <- pars2[4]
  N0     <- pars2[6]

  dummy <- missnumspec #change this if you want to include missing species. it's probably not difficult, just modify pn for LM and LS

  brtsM1 <- sort(abs(brtsM), decreasing = TRUE)
  if (is.null(unlist(brtsS))) {brtsS1 <- td} else
  {
    brtsS1 <- sort(abs(c(brtsS, td)), decreasing = TRUE)
  }
  td <- abs(td) * sign(brtsM1[1])
  tsplit <- abs(tsplit) * sign(brtsM1[1])

  testit::assert(all(sign(brtsM1) == sign(td)))
  testit::assert(
    all(sign(brtsS1 * !is.null(brtsS1)) == sign(td * !is.null(brtsS1)))
  )
  testit::assert(all(sign(tsplit) == sign(td)))
  if (!(tsplit %in% brtsM1)) {stop('tsplit has to be in the main clade branching times!')}

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

  likM_pre_shift  <- sls::combine_pns_nodivision(lambda = lambdas[1], mu = mus[1], ts = tsM_pre_shift, tbar = td, nmax = nmax); log(likM_pre_shift)
  likM_post_shift <- prod(
    sls::pn(n = 1, lambda = lambdas[1], mu = mus[1], t = tsM_post_shift)
  ) * sls:::pn(n = 1, t = td, lambda = lambdas[1], mu = mus[1])^(length(tsM_pre_shift) - 1); log(likM_post_shift)
  likS_post_shift <- prod(
    sls::pn(n = 1, lambda = lambdas[2], mu = mus[2], t = brtsS1        )
  ); log(likS_post_shift)
  loglikM0 <- log(likM_pre_shift) + log(likM_post_shift)
  loglikS0 <- log(likS_post_shift)

  logcombinatoricsM <- ifelse(length(brtsM) > 1, sum(log(kvecM_before1[-1])), 0)
  logcombinatoricsS <- ifelse(length(brtsS) > 0, lfactorial(length(brtsS))   , 0)
  lM <- length(brtsM1[brtsM1 != brtsM1[1]]) # number of speciations in the Main clade
  lS <- length(brtsS1[brtsS1 != brtsS1[1]]) # number of speciations in the Sub clade
  loglikM <- loglikM0 + logcombinatoricsM + log(lambdas[1] + (length(brtsM) == 0)) * lM
  loglikS <- loglikS0 + logcombinatoricsS + log(lambdas[2] + (length(brtsS) == 0)) * lS

  Pc <- 1
  if (!missing(cond))
  {
    if (cond %in% c(1,2,3))
    {
      conditioning <- sls::Pc_1shift(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS)
      Pc <- conditioning[[cond]]; Pc
    }
  }

  loglik <- loglikM + loglikS - log(Pc); loglik
  return(loglik)
}

#' @title Q-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood integrating the Q-equation
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_slsQ <- function(pars1, pars2 = c(100, 1, 1, brtsM[2], 0, 2), brtsM, brtsS)
{
  missnumspec <- c(0,0)
  lambdas <- c(pars1[1], pars1[4])
  mus     <- c(pars1[2], pars1[5])
  Ks      <- c(pars1[3], pars1[6])
  td      <- pars1[7]

  lx     <- pars2[1]
  cond   <- pars2[3]
  tsplit <- pars2[4]
  N0s    <- c(pars2[6], 1)

  testit::assert(all(sign(brtsM) == sign(td)))
  testit::assert(
    all(sign(brtsS * !is.null(brtsS)) == sign(td * !is.null(brtsS)))
  )
  testit::assert(all(sign(tsplit) == sign(td)))
  if (!(tsplit %in% brtsM)) {stop('tsplit has to be in the main clade branching times!')}

  #BASIC SETTINGS AND CHECKS
  Nclades <- length(lambdas)
  brtsM1 <- sort(c(0, abs(c(brtsM, td))), decreasing = TRUE)
  brtsS1 <- sort(c(0, abs(c(brtsS, td))), decreasing = TRUE)
  brts_list <- list(brtsM = brtsM1, brtsS = brtsS1)
  abstol <- 1e-16; reltol <- 1e-10
  nvec <- 0:lx
  logliks <- rep(NA, Nclades)

  #LIKELIHOOD INTEGRATION
  clade <- 0 #clade == 1 is the main clade, clade == 2 is the subclade
  while ((clade <- clade + 1) <= Nclades)
  {
    #SETTING CLADE CONDITIONS
    lambda <- lambdas[clade]
    mu     <- mus[clade]
    K      <- Ks[clade]
    N0     <- N0s[clade]
    maxT   <- length(brts_list[[clade]])
    brts   <- brts_list[[clade]]

    #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
    Qi <- c(1, rep(0, lx))
    Qt <- matrix(0, ncol = (lx + 1), nrow = maxT)
    Qt[1,] <- Qi
    dimnames(Qt)[[2]] <- paste0("Q", 0:lx)
    k <- N0
    t <- 2
    D <- C <- rep(1, maxT)

    #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
    while (t <= maxT)
    {
      # print(t)
      # print(abs(brts[t] - brts[t - 1]))
      #Applying A operator
      if (lambda == 0 && mu == 0)
      {
        Qt[t,] <- Qt[(t-1),]
      }else
      {
        transition_matrix <- DDD:::dd_loglik_M_aux(pars = c(lambda, mu, K), lx = lx + 1, k = k, ddep = 1)
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
          Qt[t,] <- Qt[t,] * k * lambda
          k <- k + 1
        }else
        {
          # print("-1")
          Qt[t,] <- Qt[t,] * (k + nvec)^-1
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

  Pc <- 1
  if (!missing(cond))
  {
    if (cond %in% c(1,2,3))
    {
      conditioning <- sls::Pc_1shift(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS)
      Pc <- conditioning[[cond]]; Pc
    }
  }

  total_loglik <- sum(logliks) - log(Pc)
  return(total_loglik)
}

#' @title DDD-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood calling the routine from the DDD package
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_DDD <- function(pars1, pars2 = c(100, 1, 1, brtsM[2], 0, 2), brtsM, brtsS, missnumspec = c(0,0))
{
  # lambdas <- c(pars1[1], pars1[4])
  # mus     <- c(pars1[2], pars1[5])
  # Ks      <- c(pars1[3], pars1[6])
  # td      <- pars1[7]
  # nmax   <- pars2[1]

  # tsplit <- pars2[4]
  # N0     <- pars2[6]
  cond <- pars2[3]
  # pars2[3] <- 0 # i will impose my conditioning
  pars1copy <- pars1; pars1copy[7] <- abs(pars1copy[7])
  pars2copy <- pars2; pars2copy[4] <- abs(pars2copy[4]); pars2copy[3] <- 0 #i will impose my conditioning
  testit::assert(pars2copy[4] %in% abs(brtsM)) #tsplit in brtsM
  loglik <- DDD::dd_KI_loglik(pars1 = pars1copy, pars2 = pars2copy,
                              brtsM = abs(brtsM), brtsS = abs(brtsS), missnumspec = missnumspec)

  Pc <- 1
  if (!missing(cond))
  {
    if (cond %in% c(1,2,3))
    {
      pars2[3] <- cond
      conditioning <- sls::Pc_1shift(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS)
      Pc <- conditioning[[cond]]; Pc
    }
  }

  loglikC <- loglik - log(Pc)
  return(loglikC)
}
