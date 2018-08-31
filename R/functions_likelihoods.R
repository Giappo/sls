#' @title P-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convoluting Nee's functions
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
lik_shift_P <- function(pars, brtsM, brtsS, nmax = 1e2, cond = 0) {

  lambdas <- c(pars[1], pars[3])
  mus     <- c(pars[2], pars[4])

  tshift <- brtsS[1]
  if (length(tshift) > 0)
  {
    vec1 <- matrix(rbind(brtsM[brtsM < tshift], rep(1, length(brtsM[brtsM < tshift]))), nrow = 2)
    vec2 <- matrix(rbind(tshift, rep(-1, length(tshift))), nrow = 2)
    vec3 <- matrix(rbind(brtsM[brtsM > tshift], rep(1, length(brtsM[brtsM > tshift]))), nrow = 2)

    M  <- cbind(vec1, vec2, vec3)
    M2 <- M
    M2[2,] <- (cumsum(M[2,]) - 1)
    M3 <- M2[,-sign(M[2,]) < 0]
    dim(M3) <- c(2, sum(sign(M[2,]) >= 0))
    kvec <- M3[2,]

    tsM_preshift  <- tshift - brtsM[brtsM < tshift]
    tsM_postshift <- 0 - brtsM[brtsM > tshift]; if (length(tsM_postshift) == 0) {tsM_postshift <- 0}
    tbar <- 0 - tshift

    if (length(tsM_preshift) == 0) {cat("There are no branching times before the shift"); return(-Inf)}

    likM_preshift  <- sls::combine_pns(lambda = lambdas[1], mu = mus[1], ts = tsM_preshift, tbar = tbar, nmax = nmax)
    likM_postshift <- prod(
      sls::pn(n = 1, t = tsM_postshift, lambda = lambdas[1], mu = mus[1])
    ) * sls:::pn(n = 1, t = tbar, lambda = lambdas[1], mu = mus[1])^(length(tsM_preshift) - 1)
    likS_postshift <- prod(
      sls::pn(n = 1, t = abs(brtsS), lambda = lambdas[2], mu = mus[2])
    )
    loglikM0 <- log(likM_preshift) + log(likM_postshift)
    loglikS0 <- log(likS_postshift)
  }else
  {
    M <- rbind(brtsM, rep(1, length(brtsM)))
    M2 <- M
    M2[2,] <- (cumsum(M[2,]) - 1)
    M3 <- M2[,-sign(M[2,]) < 0]
    dim(M3) <- c(2, sum(sign(M[2,]) >= 0))
    kvec <- M3[2,]

    likM <- prod(sls::pn(n = 1, t = 0 - brtsM, lambda = lambdas[1], mu = mus[1]))
    loglikM0 <- log(likM)
    loglikS0 <- 0
  }

  logcombinatoricsM <- ifelse(length(kvec ) > 1, sum(log(kvec[2:length(kvec)])), 0)
  logcombinatoricsS <- ifelse(length(brtsS) > 0, lfactorial(length(brtsS) - 1) , 0)
  loglikM <- loglikM0 + (length(brtsM[brtsM != brtsM[1]])) * log(lambdas[1]) + logcombinatoricsM
  loglikS <- loglikS0 + logcombinatoricsS
  if ((length(brtsS) - 1) > 0) {loglikS <- loglikS + (length(brtsS) - 1) * log(lambdas[2])} #adding lambda terms

  Pc <- 1
  if (!missing(cond))
  {
    if (cond %in% c(1,2,3))
    {
      conditioning <- sls::Pc_1shift(pars = pars, brtsM = brtsM, brtsS = brtsS)
      Pc <- conditioning[[cond]]; Pc
    }
  }

  # loglik <- log(likM) + log(likS) - log(Pc); loglik
  loglik <- loglikM + loglikS - log(Pc); loglik

  return(loglik)
}

#' @title P-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convoluting Nee's functions
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
lik_shift_P2 <- function(pars1, pars2 = c(100, 1, 1, brtsM[2], 0, 2), brtsM, brtsS, missnumspec = c(0,0)) {

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
  brtsS1 <- sort(abs(c(brtsS, td)), decreasing = TRUE)
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
      conditioning <- sls::Pc_1shift2(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS)
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
lik_shift_Q <- function(pars,
                        brtsM,
                        brtsS,
                        lx = floor(min(20 * max(length(brtsM), length(brtsS)), 1000)),
                        cond = 0
){

  lambdas <- c(pars[1], pars[3])
  mus     <- c(pars[2], pars[4])
  missnumspec = c(0,0)

  #BASIC SETTINGS AND CHECKS
  Nclades <- length(lambdas)
  brts_list <- list(brtsM = brtsM, brtsS = brtsS)
  shift_times <- unlist(lapply(brts_list, FUN = function(x) x[1]))
  abstol <- 1e-16; reltol <- 1e-10

  #ADJUSTING DATA
  nvec <- 0:lx
  clade <- 0 #clade == 1 is the main clade, clade == 2 is the subclade
  logliks <- rep(NA, Nclades)
  #LIKELIHOOD INTEGRATION
  while ((clade <- clade + 1) <= Nclades)
  {
    #SETTING CLADE CONDITIONS
    shift_times2 <- shift_times[shift_times > min(brts_list[[clade]])]
    time_points <- sort(unique(c(brts_list[[clade]], shift_times2)), decreasing = FALSE)
    # data <- MBD:::brts2time_intervals_and_births(time_points); time_intervals <- data$time_intervals; time_intervals <- c(0, time_intervals)
    time_intervals <- sls::brts2time_intervals(time_points)

    lambda <- lambdas[clade]
    mu     <- mus[clade]
    N0     <- sum(brts_list[[clade]] == brts_list[[clade]][1])

    #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
    Qi <- c(1, rep(0, lx))
    Qt <- matrix(0, ncol = (lx + 1), nrow = (length(time_intervals)))
    Qt[1,] <- Qi
    dimnames(Qt)[[2]] <- paste0("Q", 0:lx)
    k <- N0
    t <- 2
    D <- C <- rep(1, (length(time_intervals)))

    #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
    while (t <= length(time_intervals))
    {
      #Applying A operator
      if (lambda == 0 && mu == 0)
      {
        Qt[t,] <- Qt[(t-1),]
      }else
      {
        transition_matrix <- DDD:::dd_loglik_M_aux(pars = c(lambda, mu, Inf), lx = lx + 1, k = k, ddep = 1)
        Qt[t,] <- abs(expoRkit::expv(v = Qt[(t-1),], x = transition_matrix, t = time_intervals[t]))
      }
      # if (methode == "analytical") {
      # }
      # else {
      #   transition_matrix <- MBD:::create_A(lambda = lambda, mu = mu, nu = 0, q = 0, k = k, max_number_of_species = lx)
      #   Qt[t,] <- deSolve::ode(y = Qt[(t-1),], times = c(0, time_intervals[t]), func = MBD:::mbd_loglik_rhs,
      #                          parms = transition_matrix, atol = 1e-16, rtol = 1e-10)[2,-1]
      # }

      #Applying C operator (this is a trick to avoid precision issues)
      C[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * C[t]

      #what time is it?
      tempo <- time_points[t]

      if (t < length(time_intervals))
      {
        #Applying B operator
        if (all(tempo != shift_times))
        {
          Qt[t,] <- Qt[t,] * k * lambda
          k <- k + 1
        }else
        {
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
      conditioning <- sls::Pc_1shift(pars = pars, brtsM = brtsM, brtsS = brtsS)
      Pc <- conditioning[[cond]]; Pc
    }
  }

  total_loglik <- sum(logliks) - log(Pc)
  return(total_loglik)
}

#' @title Q-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood integrating the Q-equation
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
lik_shift_Q2 <- function(pars1, pars2 = c(100, 1, 1, brtsM[2], 0, 2), brtsM, brtsS)
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
      conditioning <- sls::Pc_1shift2(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS)
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
lik_shift_DDD <- function(pars, brtsM, brtsS, cond = 0,
                          lx = floor(min(20 * max(length(brtsM), length(brtsS)), 1000)), methode = "lsoda")
{
  missnumspec = c(0,0)
  pars1 <- c(pars[1], pars[2], Inf, pars[3], pars[4], Inf, abs(brtsS[1]))
  pars2 <- rep(0, 6)
  pars2[1] <- lx
  pars2[2] <- 1 #ddep
  pars2[3] <- 0
  pars2[4] <- abs(brtsS[1])
  pars2[5] <- 0
  pars2[6] <- sum(brtsM == brtsM[1])
  brtsS    <- brtsS[-1]

  if (length(pars2) == 4) {
    pars2[5] = 0
    pars2[6] = 2
  }
  abstol = 1e-16
  reltol = 1e-14
  m = missnumspec
  brts = -sort(abs(c(brtsM, brtsS)), decreasing = TRUE)
  if (sum(brts == 0) == 0) {
    brts[length(brts) + 1] = 0
  }
  soc = pars2[6]
  S = length(brts) + (soc - 2)
  brtsM = -sort(abs(brtsM), decreasing = TRUE)
  if (sum(brtsM == 0) == 0) {
    brtsM[length(brtsM) + 1] = 0
  }
  brtsS = -sort(abs(brtsS), decreasing = TRUE)
  if (sum(brtsS == 0) == 0) {
    brtsS[length(brtsS) + 1] = 0
  }
  if (min(pars1) < 0 | -pars1[7] <= min(brtsM) | -pars1[7] >=
      min(brtsS)) {
    loglik = -Inf
  }
  else {
    if (((pars1[2] == 0 || pars1[4] == 0) && pars2[2] ==
         2) | ((pars1[1] == 0 | pars1[3] == 0) & pars2[2] ==
               4) | pars1[1] <= pars1[2] | pars1[4] <= pars1[5]) {
      cat("These parameter values cannot satisfy lambda(N) = mu(N) for some finite N.\n")
      loglik = -Inf
    }
    else {
      laM = pars1[1]
      muM = pars1[2]
      KM = pars1[3]
      laS = pars1[4]
      muS = pars1[5]
      KS = pars1[6]
      tinn = -pars1[7]
      lmax = pars2[1]
      ddep = pars2[2]
      if (ddep == 1) {
        lxM = min(max(1 + m[1], 1 + ceiling(laM/(laM -
                                                   muM) * KM)), ceiling(lmax))
        lxS = min(max(1 + m[1], 1 + ceiling(laS/(laS -
                                                   muS) * KS)), ceiling(lmax))
      }
      else if (ddep == 1.3) {
        lxM = min(max(1 + m[1], 1 + ceiling(KM)), ceiling(lmax))
        lxS = min(max(1 + m[1], 1 + ceiling(KS)), ceiling(lmax))
      }
      else {
        lxM = round(lmax)
        lxS = round(lmax)
      }
      n0 = (ddep == 2 | ddep == 4)
      # cond = pars2[3]
      tsplit = -pars2[4]
      S1 = length(brtsM) - 1 + (soc - 2)
      if (sum(brtsS == tinn) == 0) {
        brtsS = c(tinn, brtsS)
      }
      S2 = length(brtsS) - 1
      S1a = S1
      S2a = S2
      summ = sum(m)
      if (length(m) == 2) {
        S1a = S1 + m[1]
        S2a = S2 + m[2]
        summ = 0
      }
      if ((ddep == 1 & ((ceiling(laM/(laM - muM) * KM) <
                         S1a) | (ceiling(laS/(laS - muS) * KS) < S2a))) |
          (ddep == 1.3 & ((ceiling(KM) < S1a) | (ceiling(KS) <
                                                 S2a) | (ceiling(KM) + ceiling(KS) < S1a +
                                                         S2a + summ)))) {
        loglik = -Inf
      }
      else {
        if (sum(abs(brtsM - tinn) < 1e-14) == 1) {
          tinn = tinn - 1e-08
        }
        loglikM = 0
        lx = lxM
        probs = rep(0, lx)
        probs[1] = 1
        ka = sum(brtsM < tinn)
        for (k in 2:(ka + 1)) {
          k1 = k + (soc - 2)
          t1 = brtsM[k - 1]
          t2 = min(c(tinn, brtsM[k]))
          probs = DDD::dd_loglik_M(pars1[1:3], lx, k1, ddep,
                                   tt = abs(t2 - t1), probs)
          if (t2 < tinn) {
            probs = DDD::lambdamu(0:(lx - 1) + k1, c(pars1[1:3],
                                                     0), ddep)[[1]] * probs
            sumprobs = sum(probs)
            if (sumprobs <= 0) {
              loglik = -Inf
              break
            }
            else {
              loglikM = loglikM + log(sumprobs)
            }
            probs = probs/sumprobs
          }
        }
        for (k in (ka + 1):max(ka + 1, S1 + 1)) {
          k1 = k + (soc - 2)
          t1 = max(tinn, brtsM[k - 1])
          t2 = brtsM[k]
          probs = DDD::dd_loglik_M(pars1[1:3], lx, k1 - 1,
                                   ddep, tt = abs(t2 - t1), probs)
          if (k < (S1 + 1)) {
            probs = DDD::lambdamu(0:(lx - 1) + k1 - 1, c(pars1[1:3],
                                                         0), ddep)[[1]] * probs
            sumprobs = sum(probs)
            if (sumprobs <= 0) {
              loglik = -Inf
            }
            else {
              loglikM = loglikM + log(sumprobs)
            }
            probs = probs/sumprobs
          }
        }
        if (length(m) == 1) {
          loglikM = loglikM + log(probs[1 + (0:m)])
        }
        else {
          loglikM = loglikM + log(probs[1 + m[1]])
        }
        loglikS = 0
        lx = lxS
        probs = rep(0, lx)
        probs[1] = 1
        for (k in 1:S2) {
          t1 = brtsS[k]
          t2 = brtsS[k + 1]
          probs = DDD::dd_loglik_M(pars1[4:6], lx, k, ddep,
                                   tt = abs(t2 - t1), probs)
          if (k < S2) {
            probs = DDD::lambdamu(0:(lx - 1) + k, c(pars1[4:6],
                                                    0), ddep)[[1]] * probs
            sumprobs = sum(probs)
            if (sumprobs <= 0) {
              loglik = -Inf
            }
            else {
              loglikS = loglikS + log(sumprobs)
            }
            probs = probs/sumprobs
          }
        }
        if (length(m) == 1) {
          loglikS = loglikS + log(probs[1 + (0:m)])
        }
        else {
          loglikS = loglikS + log(probs[1 + m[2]])
        }
        if (length(m) == 1) {
          loglikMmax = max(loglikM)
          loglikSmax = max(loglikS)
          loglikMdelta = loglikM - loglikMmax
          loglikSdelta = loglikS - loglikSmax
          loglik = loglikMmax + loglikSmax + log(sum(exp(loglikMdelta +
                                                           loglikSdelta[length(loglikSdelta):1])))
        }
        else {
          loglik = loglikM + loglikS
        }
        if (is.nan(loglik) | is.na(loglik)) {
          loglik = -Inf
        }
        if (cond == 0 | loglik == -Inf) {
          logliknorm = 0
        }
        else {
          tcrown = brts[1]
          tpres = 0
          lx = lxS
          probs = rep(0, lx)
          probs[2] = 1
          probs = DDD::dd_loglik_M(pars1[4:6], lx, 0, ddep,
                                   tt = abs(tpres - tinn), probs)
          PS = 1 - probs[1]
          lx = lxM
          probs = matrix(0, lx, lx)
          probs[2, 2] = 1
          dim(probs) = c(lx * lx, 1)
          ly = lx^2
          probs = DDD::dd_loglik_M2(pars = pars1[1:3], lx = lx,
                                    ddep = ddep, tt = abs(tinn - tcrown), p = probs)
          dim(probs) = c(lx, lx)
          probs[1, 1:lx] = 0
          probs[1:lx, 1] = 0
          auxM1 = rep(0:(lx - 1), times = lx) + rep(0:(lx -
                                                         1), each = lx)
          probs = probs * rep(0:(lx - 1), times = lx)/auxM1
          dim(probs) = c(lx, lx)
          probs = rbind(probs[2:lx, 1:lx], rep(0, lx))
          dim(probs) = c(lx * lx, 1)
          probs = DDD::dd_loglik_M2(pars = pars1[1:3], lx = lx,
                                    ddep = ddep, tt = abs(tpres - tinn), p = probs)
          dim(probs) = c(lx, lx)
          PM12 = sum(probs[2:lx, 2:lx])
          PM2 = sum(probs[1, 2:lx])
          logliknorm = log(2) + log(PM12 + PS * PM2)
        }
        if (length(m) > 1) {
          Sv = c(S1, S2)
        }
        else {
          Sv = S
        }
        loglik = loglik - logliknorm - sum(lgamma(Sv +
                                                    m + 1) - lgamma(Sv + 1) - lgamma(m + 1))
      }
    }
  }
  if (pars2[5] == 1) {
    s1 = sprintf("Parameters: %f %f %f %f %f %f %f, ", pars1[1],
                 pars1[2], pars1[3], pars1[4], pars1[5], pars1[6],
                 pars1[7])
    s2 = sprintf("Loglikelihood: %f", loglik)
    cat(s1, s2, "\n", sep = "")
    flush.console()
  }
  loglik = as.numeric(loglik)
  if (is.nan(loglik) | is.na(loglik)) {
    loglik = -Inf
  }

  Pc <- 1
  if (!missing(cond))
  {
    if (cond %in% c(1,2,3))
    {
      conditioning <- sls::Pc_1shift(pars = pars, brtsM = brtsM, brtsS = brtsS)
      Pc <- conditioning[[cond]]; Pc
    }
  }
  loglikC <- loglik - log(Pc)
  return(loglikC)
}

#' @title DDD-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood calling the routine from the DDD package
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
lik_shift_DDD2 <- function(pars1, pars2 = c(100, 1, 1, brtsM[2], 0, 2), brtsM, brtsS, missnumspec = c(0,0))
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
  pars2copy <- pars2; pars2copy[3] <- 0 #i will impose my conditioning

  loglik <- DDD::dd_KI_loglik(pars1 = pars1, pars2 = pars2copy, brtsM = brtsM, brtsS = brtsS, missnumspec = missnumspec)

  Pc <- 1
  if (!missing(cond))
  {
    if (cond %in% c(1,2,3))
    {
      pars2[3] <- cond
      conditioning <- sls::Pc_1shift2(pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS)
      Pc <- conditioning[[cond]]; Pc
    }
  }

  loglikC <- loglik - log(Pc)
  return(loglikC)
}

#' @title DDD-conditioning
#' @author Giovanni Laudanno
#' @description Calculates the conditioning from the DDD package for the Key Innovation model
#' @inheritParams default_params_doc
#' @return Conditioning probability
#' @export
DDD_conditioning <- function(pars, brtsM, brtsS, lx, ddep = 1) {

  lambdas <- c(pars[1], pars[3])
  mus     <- c(pars[2], pars[4])
  pars1 <- c(lambdas[1], mus[1], Inf, lambdas[2], mus[2], Inf, brtsS[1])
  laM = pars1[1]
  muM = pars1[2]
  KM = pars1[3]
  laS = pars1[4]
  muS = pars1[5]
  KS = pars1[6]
  tinn = -abs(pars1[7])
  lmax = lx #pars2[1]
  # ddep = pars2[2]

  m = 0
  lxM = min(max(1 + m[1], 1 + ceiling(KM)), ceiling(lmax))
  lxS = min(max(1 + m[1], 1 + ceiling(KS)), ceiling(lmax))

  tcrown = brtsM[1]
  tpres  = 0
  # tinn <- brtsS[1]
  lx = lxS
  # lxS <- lxM <- lx
  probs = rep(0, lx)
  probs[2] = 1
  probs = DDD:::dd_loglik_M(pars1[4:6], lx, 0, ddep,
                            tt = abs(tpres - tinn), probs)
  PS = 1 - probs[1]
  lx = lxM
  probs = matrix(0, lx, lx)
  probs[2, 2] = 1
  dim(probs) = c(lx * lx, 1)
  ly = lx^2
  probs = DDD:::dd_loglik_M2(pars = pars1[1:3], lx = lx,
                             ddep = ddep, tt = abs(tinn - tcrown), p = probs)
  dim(probs) = c(lx, lx)
  probs[1, 1:lx] = 0
  probs[1:lx, 1] = 0
  auxM1 = rep(0:(lx - 1), times = lx) + rep(0:(lx - 1), each = lx)
  probs = probs * rep(0:(lx - 1), times = lx)/auxM1
  dim(probs) = c(lx, lx)
  probs = rbind(probs[2:lx, 1:lx], rep(0, lx))
  dim(probs) = c(lx * lx, 1)
  probs = DDD:::dd_loglik_M2(pars = pars1[1:3], lx = lx,
                             ddep = ddep, tt = abs(tpres - tinn), p = probs)
  dim(probs) = c(lx, lx)
  PM12 = sum(probs[2:lx, 2:lx])
  PM2 = sum(probs[1, 2:lx])
  logliknorm = log(2) + log(PM12 + PS * PM2); exp(logliknorm)
  return(exp(logliknorm))
}

#' @title sls-conditioning
#' @author Giovanni Laudanno
#' @description Calculates three different kind of conditioning probabilities
#' @inheritParams default_params_doc
#' @return Three different conditioning probabilities
#' @export
Pc_1shift <- function(pars, brtsM, brtsS, Nmax = 100) {

  lambdas <- c(pars[1], pars[3])
  mus     <- c(pars[2], pars[4])

  if (length(brtsS) > 0)
  {
    tc <- brtsM[1]
    ts <- brtsS[1]
    tp <- 0
    A <- abs(ts - tc); B <- abs(tp - ts)
    PS <- 1 - pn(n = 0, t = B, lambda = lambdas[2], mu = mus[2])
    nvec <- 1:Nmax
    ns1  <- row(matrix(NA, nrow = Nmax, ncol = Nmax))
    ns2  <- col(matrix(NA, nrow = Nmax, ncol = Nmax))
    pA   <- sls::pt(t = A, lambda = lambdas[1], mu = mus[1]); pA
    uA   <- sls::ut(t = A, lambda = lambdas[1], mu = mus[1]); uA
    pB1  <- sls::pt(t = B, lambda = lambdas[1], mu = mus[1]); pB1
    pB2  <- sls::pt(t = B, lambda = lambdas[2], mu = mus[2]); pB2
    pns1 <- sls::pn(n = ns1, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns1) <- paste0("ns1=", nvec); colnames(pns1) <- paste0("ns2=", nvec); head(pns1)
    pns2 <- sls::pn(n = ns2, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns2) <- paste0("ns1=", nvec); colnames(pns2) <- paste0("ns2=", nvec); head(pns2)
    aux1 <- pns1 * pns2 * (ns1/(ns1 + ns2)) * (1 - (1 - pB1)^ns2)
    P1   <- sum(aux1) #branch 2 survives till the present
    aux2 <- aux1 * (1 - (1 - pB1)^(ns1 - 1)); head(aux2)
    P2   <- sum(aux2) #both branches 1 and 2 survive till the present
  }else
  {
    P1 <- (1 - sls::pn(n = 0, lambda = lambdas[1], mu = mus[1]))
    P2 <- (1/2) * (1 - sls::pn(n = 0, lambda = lambdas[1], mu = mus[1]))^2 #both branches 1 and 2 survive till the present
    PS <- 0
  }

  Pc1  <- 2 * PS * P1 + 2 * (1 - PS) * P2
  Pc2  <- 2 * PS * P2
  Pc3  <- 2 * PS * P1

  return(list(Pc1 = Pc1, Pc2 = Pc2, Pc3 = Pc3))
}

#' @title sls-conditioning
#' @author Giovanni Laudanno
#' @description Calculates three different kind of conditioning probabilities
#' @inheritParams default_params_doc
#' @return Three different conditioning probabilities
#' @export
Pc_1shift2 <- function(pars1, pars2, brtsM, brtsS) {

  missnumspec <- c(0,0)
  lambdas <- c(pars1[1], pars1[4])
  mus     <- c(pars1[2], pars1[5])
  Ks      <- c(pars1[3], pars1[6])
  td      <- pars1[7]
  td      <- abs(td) * sign(brtsM[1])

  lx     <- pars2[1]
  cond   <- pars2[3]
  tsplit <- pars2[4]
  soc    <- pars2[6]

  if (soc != 2) {stop("Pc can be calculated only if phylogeny starts with a crown!")}

  if (length(brtsS) > 0)
  {
    tc <- brtsM[1]
    ts <- td
    tp <- 0
    A <- abs(ts - tc); B <- abs(tp - ts)
    PS <- 1 - pn(n = 0, t = B, lambda = lambdas[2], mu = mus[2])
    nvec <- 1:lx
    ns1  <- row(matrix(NA, nrow = lx, ncol = lx))
    ns2  <- col(matrix(NA, nrow = lx, ncol = lx))
    pA   <- sls::pt(t = A, lambda = lambdas[1], mu = mus[1]); pA
    uA   <- sls::ut(t = A, lambda = lambdas[1], mu = mus[1]); uA
    pB1  <- sls::pt(t = B, lambda = lambdas[1], mu = mus[1]); pB1
    pB2  <- sls::pt(t = B, lambda = lambdas[2], mu = mus[2]); pB2
    pns1 <- sls::pn(n = ns1, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns1) <- paste0("ns1=", nvec); colnames(pns1) <- paste0("ns2=", nvec); head(pns1)
    pns2 <- sls::pn(n = ns2, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns2) <- paste0("ns1=", nvec); colnames(pns2) <- paste0("ns2=", nvec); head(pns2)
    aux1 <- pns1 * pns2 * (ns1/(ns1 + ns2)) * (1 - (1 - pB1)^ns2)
    P1   <- sum(aux1) #branch 2 survives till the present
    aux2 <- aux1 * (1 - (1 - pB1)^(ns1 - 1)); head(aux2)
    P2   <- sum(aux2) #both branches 1 and 2 survive till the present
  }else
  {
    P1 <- (1 - sls::pn(n = 0, lambda = lambdas[1], mu = mus[1]))
    P2 <- (1/2) * (1 - sls::pn(n = 0, lambda = lambdas[1], mu = mus[1]))^2 #both branches 1 and 2 survive till the present
    PS <- 0
  }

  Pc1  <- 2 * PS * P1 + 2 * (1 - PS) * P2
  Pc2  <- 2 * PS * P2
  Pc3  <- 2 * PS * P1

  return(list(Pc1 = Pc1, Pc2 = Pc2, Pc3 = Pc3))
}

#' #OLD SCRIPTS
#' #' @title P-likelihood
#' #' @author Giovanni Laudanno
#' #' @description Calculates the likelihood convoluting Nee's functions
#' #' @inheritParams default_params_doc
#' #' @return The likelihood
#' #' @export
#' lik_shift_P2 <- function(pars, LM, LS, nmax = 1e2, cond = 0) {
#'
#'   lambdas <- c(pars[1], pars[3])
#'   mus     <- c(pars[2], pars[4])
#'
#'   tshift <- LS[1,1]
#'   LTTS <- sls::L2LTT(LS, shifts = NULL)$LTT_reconstructed;   LS; LTTS
#'   LTTM <- sls::L2LTT(LM, shifts = tshift)$LTT_reconstructed; LM; LTTM
#'   kvecM <- LTTM[2, c(1, diff(LTTM[2,])) > 0]
#'   brtsM <- LTTM[1, c(1, diff(LTTM[2,])) > 0]
#'   if (ncol(LTTS) > 0)
#'   {
#'     brtsS <- LTTS[1, c(1, diff(LTTS[2,])) > 0]
#'   }else
#'   {
#'     brtsS <- 0
#'   }
#'   if (brtsS[1] > 0)
#'   {
#'     tsM_preshift  <- brtsM[brtsM > tshift] - tshift
#'     tsM_postshift <- brtsM[brtsM < tshift]
#'     tbar <- tshift
#'     if (length(tsM_postshift) == 0) {tsM_postshift <- 0}
#'     if (length(tsM_preshift) == 0) {cat("There are no branching times before the shift"); return(-Inf)}
#'
#'     likM_preshift  <- sls::combine_pns(lambda = lambdas[1], mu = mus[1], ts = tsM_preshift, tbar = tbar, nmax = nmax)
#'     likM_postshift <- prod(
#'       sls::pn(n = 1, t = tsM_postshift, lambda = lambdas[1], mu = mus[1])
#'     ) * sls:::pn(n = 1, t = tbar, lambda = lambdas[1], mu = mus[1])^(length(tsM_preshift) - 1)
#'     likS_postshift <- prod(
#'       sls::pn(n = 1, t = abs(brtsS), lambda = lambdas[2], mu = mus[2])
#'     )
#'     loglikM0 <- log(likM_preshift) + log(likM_postshift)
#'     loglikS0 <- log(likS_postshift)
#'   }else
#'   {
#'     likM <- prod(sls::pn(n = 1, t = 0 - brtsM, lambda = lambdas[1], mu = mus[1]))
#'     loglikM0 <- log(likM)
#'     loglikS0 <- 0
#'   }
#'
#'   logcombinatoricsM <- ifelse(length(kvecM) > 1, sum(log(kvecM[2:length(kvecM)])), 0)
#'   logcombinatoricsS <- ifelse(length(brtsS) > 0, lfactorial(length(brtsS) - 1)   , 0)
#'   loglikM <- loglikM0 + (length(brtsM[brtsM != brtsM[1]])) * log(lambdas[1]) + logcombinatoricsM
#'   loglikS <- loglikS0 + logcombinatoricsS
#'   if ((length(brtsS) - 1) > 0) {loglikS <- loglikS + (length(brtsS) - 1) * log(lambdas[2])} #adding lambda terms
#'
#'   Pc <- 1
#'   if (!missing(cond))
#'   {
#'     if (cond %in% c(1,2,3))
#'     {
#'       conditioning <- sls::Pc_1shift(pars = pars, brtsM = brtsM, brtsS = brtsS)
#'       Pc <- conditioning[[cond]]; Pc
#'     }
#'   }
#'
#'   loglik <- loglikM + loglikS - log(Pc); loglik
#'   return(loglik)
#' }
#'
#' #' @title Q-likelihood
#' #' @author Giovanni Laudanno
#' #' @description Calculates the likelihood integrating the Q-equation
#' #' @inheritParams default_params_doc
#' #' @return The likelihood
#' #' @export
#' lik_shift_Q2 <- function(pars,
#'                          LM,
#'                          LS,
#'                          lx = floor(min(20 * max(length(brtsM), length(brtsS)), 1000)),
#'                          cond = 0
#' ){
#'
#'   lambdas <- c(pars[1], pars[3])
#'   mus     <- c(pars[2], pars[4])
#'   missnumspec = c(0,0)
#'
#'   tshift <- LS[1,1]
#'   LTTS <- sls::L2LTT(LS, shifts = NULL)$LTT_reconstructed;   LS; LTTS
#'   LTTM <- sls::L2LTT(LM, shifts = tshift)$LTT_reconstructed; LM; LTTM
#'   LTT.List <- list(LTTM, LTTS)
#'   kvecM <- LTTM[2, c(1, diff(LTTM[2,])) > 0]
#'   brtsM <- LTTM[1, c(1, diff(LTTM[2,])) > 0]
#'   if (ncol(LTTS) > 0)
#'   {
#'     brtsS <- LTTS[1, c(1, diff(LTTS[2,])) > 0]
#'   }else
#'   {
#'     brtsS <- 0
#'   }
#'
#'   #BASIC SETTINGS AND CHECKS
#'   Nclades <- length(lambdas)
#'   brts_list <- list(brtsM = brtsM, brtsS = brtsS)
#'   shift_times <- unlist(lapply(brts_list, FUN = function(x) x[1]))
#'   abstol <- 1e-16; reltol <- 1e-10
#'
#'   #ADJUSTING DATA
#'   nvec <- 0:lx
#'   clade <- 0 #clade == 1 is the main clade, clade == 2 is the subclade
#'   logliks <- rep(NA, Nclades)
#'   #LIKELIHOOD INTEGRATION
#'   while ((clade <- clade + 1) <= Nclades)
#'   {
#'     #SETTING CLADE CONDITIONS
#'     # shift_times2 <- shift_times[shift_times > min(brts_list[[clade]])]
#'     # time_points <- sort(unique(c(brts_list[[clade]], shift_times2)), decreasing = FALSE)
#'     # data <- MBD:::brts2time_intervals_and_births(time_points); time_intervals <- data$time_intervals; time_intervals <- c(0, time_intervals)
#'     time_points <- LTT.List[[clade]][1,]
#'     kvec        <- LTT.List[[clade]][2,]
#'     kdiff       <- c(1, diff(kvec))
#'     # time_intervals <- sls::brts2time_intervals(time_points)
#'
#'     lambda <- lambdas[clade]
#'     mu     <- mus[clade]
#'     N0     <- kvec[1]
#'
#'     #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
#'     Qi <- c(1, rep(0, lx))
#'     Qt <- matrix(0, ncol = (lx + 1), nrow = (length(time_points)))
#'     Qt[1,] <- Qi
#'     dimnames(Qt)[[2]] <- paste0("Q", 0:lx)
#'     k <- N0
#'     t <- 2
#'     D <- C <- rep(1, (length(time_points)))
#'
#'     #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
#'     while (t <= length(time_intervals))
#'     {
#'       #Applying A operator
#'       if (lambda == 0 && mu == 0)
#'       {
#'         Qt[t,] <- Qt[(t-1),]
#'       }else
#'       {
#'         transition_matrix <- DDD:::dd_loglik_M_aux(pars = c(lambda, mu, Inf), lx = lx + 1, k = k, ddep = 1)
#'         Qt[t,] <- abs(expoRkit::expv(v = Qt[(t-1),], x = transition_matrix, t = abs(time_points[t] - time_points[t - 1])))
#'       }
#'       # if (methode == "analytical") {
#'       # }
#'       # else {
#'       #   transition_matrix <- MBD:::create_A(lambda = lambda, mu = mu, nu = 0, q = 0, k = k, max_number_of_species = lx)
#'       #   Qt[t,] <- deSolve::ode(y = Qt[(t-1),], times = c(0, time_intervals[t]), func = MBD:::mbd_loglik_rhs,
#'       #                          parms = transition_matrix, atol = 1e-16, rtol = 1e-10)[2,-1]
#'       # }
#'
#'       #Applying C operator (this is a trick to avoid precision issues)
#'       C[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * C[t]
#'
#'       #what time is it?
#'       tempo <- time_points[t]
#'
#'       if (t < length(time_points))
#'       {
#'         #Applying B operator
#'         # if (all(tempo != shift_times))
#'         if (kdiff[t] == 1)
#'         {
#'           Qt[t,] <- Qt[t,] * k * lambda
#'           k <- k + 1
#'         }else if (kdiff[t] == -1)
#'         {
#'           Qt[t,] <- Qt[t,] * (k + nvec)^-1
#'           k <- k - 1
#'         }
#'
#'         #Applying D operator (this works exactly like C)
#'         D[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * D[t]
#'
#'         #Updating running parameters
#'         t <- t + 1
#'       }else{break}
#'     }
#'
#'     #Selecting the state I am interested in
#'     vm <- 1/choose((k + missnumspec[clade]), k)
#'     P  <- vm * Qt[t, (missnumspec[clade] + 1)] #I have to include +1 because of Q0
#'
#'     #Removing C and D effects from the LL
#'     loglik <- log(P) - sum(log(C)) - sum(log(D))
#'
#'     #Various checks
#'     loglik <- as.numeric(loglik)
#'     if (is.nan(loglik) | is.na(loglik))
#'     {
#'       loglik <- -Inf
#'     }
#'     logliks[clade] <- loglik
#'   }
#'
#'   Pc <- 1
#'   if (!missing(cond))
#'   {
#'     if (cond %in% c(1,2,3))
#'     {
#'       conditioning <- sls::Pc_1shift(pars = pars, brtsM = brtsM, brtsS = brtsS)
#'       Pc <- conditioning[[cond]]; Pc
#'     }
#'   }
#'
#'   total_loglik <- sum(logliks) - log(Pc)
#'   return(total_loglik)
#' }
