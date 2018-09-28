#' @title P-likelihood (old version)
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convoluting Nee's functions
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
lik_shift_P0 <- function(pars, brtsM, brtsS, nmax = 1e2, cond = 0) {

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
      conditioning <- sls::Pc_1shift0(pars = pars, brtsM = brtsM, brtsS = brtsS)
      Pc <- conditioning[[cond]]; Pc
    }
  }

  # loglik <- log(likM) + log(likS) - log(Pc); loglik
  loglik <- loglikM + loglikS - log(Pc); loglik

  return(loglik)
}

#' @title Q-likelihood (old version)
#' @author Giovanni Laudanno
#' @description Calculates the likelihood integrating the Q-equation
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
lik_shift_Q0 <- function(pars,
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

#' @title DDD-likelihood (old version)
#' @author Giovanni Laudanno
#' @description Calculates the likelihood calling the routine from the DDD package
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
lik_shift_DDD0 <- function(pars, brtsM, brtsS, cond = 0,
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

#' @title Internal MBD function
#' @description Internal MBD function.
#' @details This is not to be called by the user.
#' @export
sls_loglik_choosepar0 <- function(loglik_function = sls:::loglik_slsP0,
                                  trparsopt,
                                  trparsfix,
                                  idparsopt = 1:4,
                                  idparsfix = (1:4)[-idparsopt],
                                  brtsM,
                                  brtsS,
                                  cond = 1,
                                  pars.transform = 0
                                  # missnumspec = c(0,0)
){
  #This function provides a likelihood for a subset of parameters. This is built to work inside mbd_minusLL_vs_single_parameter or any optimizer like simplex, optim or subplex
  #idparsopt are the ids of the parameters you want to analyze
  #trparsopt are the values for parameters you want to analyze
  #idparsfix are the ids of the parameters you want to fix
  #trparsfix are the values for parameters you want to fix

  # namepars = c("la_M", "mu_M", "K_M", "la_S", "mu_S", "K_S", "t_d"); Npars <- length(namepars);
  namepars = c("la1", "mu1", "la2", "mu2"); Npars <- length(namepars);
  if (length(trparsopt) == Npars && missing(trparsfix)){trparsfix <- NULL}
  trpars1 = rep(0, Npars)
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0)
  {
    trpars1[idparsfix] <- trparsfix
  }
  if (min(trpars1[1:Npars]) < 0){loglik <- -Inf}else
  {
    if (pars.transform == 1)
    {
      #Rampal's transformation
      pars1 <- trpars1
      # pars1[1:(Npars - 1)] = trpars1[1:(Npars - 1)]/(1 - trpars1[1:(Npars - 1)]) #this for the 7 parameter configuration
      pars1 = trpars1/(1 - trpars1)
    }else
    {
      pars1 <- trpars1
    }
    loglik <- loglik_function(pars = pars1, brtsM = brtsM, brtsS = brtsS, cond = cond)
  }
  if (is.nan(loglik) || is.na(loglik))
  {
    cat("There are parameter values used which cause numerical problems:",trpars1,"\n")
    loglik <- -Inf
  }
  return(loglik)
}

#' @title Internal MBD function (old)
#' @description Internal MBD function.
#' @details This is not to be called by the user.
#' @export
sls_ML0 <- function(loglik_function = sls::loglik_slsP0,
                    initparsopt,
                    parsfix,
                    idparsopt = 1:4,
                    idparsfix = (1:4)[-idparsopt],
                    brtsM,
                    brtsS,
                    cond = 1,
                    optimmethod = 'simplex',
                    pars.transform = 1,
                    res = 10 * (1 + length(c(brtsM, brtsS))),
                    tol = c(1E-3, 1E-4, 1E-6),
                    maxiter = 1000 * round((1.25)^length(idparsopt)),
                    changeloglikifnoconv = FALSE,
                    verbose = TRUE
){# bracket#1
  # - tol = tolerance in optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - maxiter = the maximum number of iterations in the optimization
  # - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
  # - optimmethod = 'simplex' (current default) or 'subplex' (default of previous versions)
  if (!is.numeric(brtsM))
  {
    stop("'brtsM' must be numeric")
  }
  if (!is.numeric(brtsS))
  {
    stop("'brtsS' must be numeric")
  }
  if (length(idparsfix) == 0) {idparsfix <- NULL}
  if (missing(parsfix) && (length(idparsfix) == 0)){parsfix <- idparsfix <- NULL}

  options(warn=-1)
  namepars = c("la1", "mu1", "la2", "mu2"); Npars <- length(namepars); #if you add more parameters to your model just change this
  failpars <- rep(-1, Npars); names(failpars) <- namepars; #those are the parameters that you get if something goes sideways
  if (is.numeric(brtsM) == FALSE)
  {
    cat("The branching times of the main clade should be numeric.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  if (is.numeric(brtsS) == FALSE)
  {
    cat("The branching times of the sub clade should be numeric.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  idpars <- sort(c(idparsopt, idparsfix))
  if ( (sum(idpars == (1:Npars)) != Npars) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)) )
  {
    cat("The parameters to be optimized and/or fixed are incoherent.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(out2)
  }

  if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
  if (verbose == TRUE) {
    cat("You are optimizing",optstr,"\n")
  }
  if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
  if (verbose == TRUE) {
    cat("You are fixing",fixstr,"\n")
    cat("Optimizing the likelihood - this may take a while.","\n")
    flush.console()
  }
  if (pars.transform == 1)
  {
    #Rampal's transformation
    trparsopt = initparsopt/(1 + initparsopt)
    trparsopt[which(initparsopt == Inf)] = 1
    trparsfix = parsfix/(1 + parsfix)
    trparsfix[which(parsfix == Inf)] = 1
  }else
  {
    trparsopt  <- initparsopt
    trparsfix  <- parsfix
  }
  optimpars  <- c(tol, maxiter)
  initloglik <- sls:::sls_loglik_choosepar0(loglik_function = loglik_function,
                                            trparsopt = trparsopt, trparsfix = trparsfix,
                                            idparsopt = idparsopt, idparsfix = idparsfix,
                                            brtsM = brtsM, brtsS = brtsS, cond = cond,
                                            pars.transform = pars.transform)
  if (verbose == TRUE) {
    cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
    flush.console()
  }
  if (initloglik == -Inf)
  {# bracket#4
    cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return (invisible(out2))
  }
  if (verbose == TRUE) {
    sink(file = tempfile()) # Sink output here
  }
  out <- DDD:::optimizer(optimmethod = optimmethod, optimpars = optimpars,
                         fun = sls:::sls_loglik_choosepar0,
                         loglik_function = loglik_function,
                         trparsopt = trparsopt, trparsfix = trparsfix,
                         idparsopt = idparsopt, idparsfix = idparsfix,
                         brtsM = brtsM, brtsS = brtsS, cond = cond,
                         pars.transform = pars.transform)
  if (verbose == TRUE) {
    sink() # Give back the output
  }
  if (out$conv != 0)
  {# bracket#5
    cat("Optimization has not converged. Try again with different initial values.\n")
    out2 <- data.frame(t(failpars), loglik = -1, df = -1, conv = -1)
    return(invisible(out2))
  }
  MLtrpars <- as.numeric(unlist(out$par))
  if (pars.transform == 1)
  {
    #Rampal's transformation
    MLpars = MLtrpars/(1 - MLtrpars)
  }else
  {
    MLpars <- MLtrpars
  }
  MLpars1 <- rep(0, Npars); names(MLpars1) <- namepars
  MLpars1[idparsopt] <- MLpars
  if (length(idparsfix) != 0) {MLpars1[idparsfix] <- parsfix}
  ML <- as.numeric(unlist(out$fvalues))
  out2 <- data.frame(t(MLpars1), loglik = ML, df = length(initparsopt), conv = unlist(out$conv))

  tobeprint <- "Maximum likelihood parameter estimates:"
  for (ii in 1:Npars)
  {
    tobeprint <- paste(tobeprint,paste(names(MLpars1[ii]),":",sep = ""),MLpars1[ii])
  }
  if (verbose == TRUE) {
    s1 <- sprintf(tobeprint)
  }
  if(out2$conv != 0 & changeloglikifnoconv == T) {out2$loglik <- -Inf}
  if (verbose == TRUE) {
    s2 <- sprintf('Maximum loglikelihood: %f',ML)
    cat("\n",s1,"\n",s2,"\n\n")
  }

  invisible(out2)
}



#' #OLD SCRIPTS
#' #' @title P-likelihood
#' #' @author Giovanni Laudanno
#' #' @description Calculates the likelihood convoluting Nee's functions
#' #' @inheritParams default_params_doc
#' #' @return The likelihood
#' #' @export
#' lik_shift_P <- function(pars, LM, LS, nmax = 1e2, cond = 0) {
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
#' lik_shift_Q <- function(pars,
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
