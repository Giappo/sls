#' @title P-likelihood (with no division)
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convoluting Nee's functions.
#'  There is no division. It should yield the same likelihood as DDD.
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_sls_p_nodiv <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  n_max = 1e2
) {
  pars_m <- pars[1:2]
  pars_s <- pars[3:4]
  brts_m <- brts[[1]]
  brts_s <- brts[[2]]
  if (any(c(pars_m, pars_s) < 0)) {
    return(-Inf)
  }

  if (is.list(brts)) {
    n_min <- 2 * (0 + (n_0 - 1) + length(unlist(brts[[1]])))
  } else {
    n_min <- 2 * (0 + (n_0 - 1) + length(brts))
  }
  if (n_max < n_min) {
    n_max <- n_min
  }

  sls_check_input(
    brts_m = brts_m,
    brts_s = brts_s,
    cond = cond,
    n_0 = n_0,
    n_max = n_max
  )

  lambdas <- c(pars_m[1], pars_s[1])
  mus     <- c(pars_m[2], pars_s[2])
  if (any(is.infinite(c(lambdas, mus)))) {
    return(-Inf)
  }
  if (any(lambdas - mus < 0)) {
    return(-Inf)
  }

  brts_m1 <- sort(brts_m, decreasing = TRUE)
  brts_s1 <- sort(brts_s, decreasing = TRUE)
  t_d <- brts_s1[1]

  testit::assert(all(sign(brts_m1) == sign(brts_s1[1])))
  testit::assert(
    all(sign(brts_s1 * !is.null(brts_s1)) == sign(t_d * !is.null(brts_s1)))
  )

  brts_matrix <- rbind(
    brts_m1,
    rep(1, length(brts_m1))
  ); dim(brts_matrix) <- c(2, length(brts_m1))
  t_d_matrix <- c(t_d, -1); dim(t_d_matrix) <- c(2, 1)

  if (n_0 == 2) {
    brts_m2 <- c(brts_m1[1], brts_m1)
  } else {
    brts_m2 <- brts_m1
  }
  ts_m_pre_shift  <- brts_m2[brts_m2 > t_d] - t_d; ts_m_pre_shift
  ts_m_post_shift <- brts_m2[brts_m2 < t_d] ; ts_m_post_shift
  if (length(ts_m_post_shift) == 0) {
    ts_m_post_shift <- 0
  }
  if (length(ts_m_pre_shift) == 0) {
    cat("There are no branching times before the shift"); return(-Inf)
  }

  lik_m_pre_shift  <- sls::combine_pns_nodiv(
      lambda = lambdas[1],
      mu = mus[1],
      times = ts_m_pre_shift,
      tbar = t_d,
      n_max = n_max
    ); log(lik_m_pre_shift)
  log_lik_m_post_shift <- sum(
    log(
      sls::pn(
        n = 1,
        lambda = lambdas[1],
        mu = mus[1],
        t = ts_m_post_shift
      )
    )
  ) + (length(ts_m_pre_shift) - 1) *
    log(
      sls::pn(
        n = 1,
        t = t_d,
        lambda = lambdas[1],
        mu = mus[1]
      )
    )
  log_lik_s_post_shift <- sum(
    log(
      sls::pn(n = 1, lambda = lambdas[2], mu = mus[2], t = brts_s1)
    )
  )
  loglik_m0 <- log(lik_m_pre_shift) + log_lik_m_post_shift
  loglik_s0 <- log_lik_s_post_shift

  logcombinatorics_m <- logcombinatorics_s <- 0 # combinatorics

  # number of speciations in the Main clade
  l_m <- length(brts_m1[brts_m1 != brts_m1[1]])

  # number of speciations in the Subclade
  l_s <- length(brts_s1[brts_s1 != brts_s1[1]])

  loglik_m <- loglik_m0 +
    logcombinatorics_m + log(lambdas[1] + (length(brts_m) == 0)) * l_m
  loglik_s <- loglik_s0 +
    logcombinatorics_s + log(lambdas[2] + (length(brts_s) == 0)) * l_s

  pc <- sls::pc_1shift(
    pars_m = pars_m,
    pars_s = pars_s,
    brts_m = brts_m,
    brts_s = brts_s,
    cond = cond,
    n_max = n_max,
    n_0 = n_0
  )

  loglik <- loglik_m + loglik_s - log(pc)
  loglik <- as.numeric(unname(loglik))
  if (is.nan(loglik) | is.na(loglik)) {
    loglik <- -Inf
  }
  if (loglik == Inf) {
   stop("infinite loglik!")
  }
  return(loglik)
}

#' @title Q-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood integrating the Q-equation.
#'  There is no division term.
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_sls_q_nodiv <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  n_max = 1e2
) {
  pars_m <- pars[1:2]
  pars_s <- pars[3:4]
  brts_m <- brts[[1]]
  brts_s <- brts[[2]]
  if (any(c(pars_m, pars_s) < 0)) {
    return(-Inf)
  }

  if (is.list(brts)) {
    n_min <- 2 * (0 + (n_0 - 1) + length(unlist(brts[[1]])))
  } else {
    n_min <- 2 * (0 + (n_0 - 1) + length(brts))
  }
  if (n_max < n_min) {
    n_max <- n_min
  }

  sls_check_input(
    brts_m = brts_m,
    brts_s = brts_s,
    cond = cond,
    n_0 = n_0,
    n_max = n_max
  )

  lambdas <- c(pars_m[1], pars_s[1])
  mus     <- c(pars_m[2], pars_s[2])
  ks      <- c(Inf, Inf)
  if (any(is.infinite(c(lambdas, mus)))) {
    return(-Inf)
  }
  if (any(lambdas - mus < 0)) {
    return(-Inf)
  }

  brts_m1 <- sort(abs(brts_m), decreasing = TRUE)
  brts_s1 <- sort(abs(brts_s), decreasing = TRUE)
  t_d <- brts_s1[1]

  missnumspec <- c(0, 0)
  n_0s <- c(n_0, 1)

  testit::assert(all(sign(brts_m) == sign(t_d)))
  testit::assert(
    all(sign(brts_s * !is.null(brts_s)) == sign(t_d * !is.null(brts_s)))
  )

  #BASIC SETTINGS AND CHECKS
  n_clades <- length(lambdas)
  brts_m1 <- sort(c(0, abs(c(brts_m, t_d))), decreasing = TRUE)
  brts_s1 <- sort(c(0, abs(c(brts_s))), decreasing = TRUE)
  brts_list <- list(brts_m = brts_m1, brts_s = brts_s1)
  logliks <- rep(NA, n_clades)

  #LIKELIHOOD INTEGRATION
  clade <- 0 #clade == 1 is the main clade, clade == 2 is the subclade
  while (
    (clade <- clade + 1) <= n_clades
  ) {
    #SETTING CLADE CONDITIONS
    lambda <- lambdas[clade]
    mu     <- mus[clade]
    kappa  <- ks[clade]
    soc    <- n_0s[clade]
    max_t  <- length(brts_list[[clade]])
    brts   <- brts_list[[clade]]

    #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
    q_i <- c(1, rep(0, n_max))
    q_t <- matrix(0, ncol = (n_max + 1), nrow = max_t)
    q_t[1, ] <- q_i
    dimnames(q_t)[[2]] <- paste0("Q", 0:n_max)
    k <- soc
    t <- 2
    D <- C <- rep(1, max_t)

    #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
    while (t <= max_t) {
      #Applying A operator
      if (lambda == 0 && mu == 0) {
        q_t[t, ] <- q_t[(t - 1), ]
      } else {
        transition_matrix <- DDD:::dd_loglik_M_aux(
          pars = c(lambda, mu, kappa),
          lx = n_max + 1,
          k = k,
          ddep = 1
        )
        q_t[t, ] <- abs(expoRkit::expv(
          v = q_t[(t - 1), ],
          x = transition_matrix,
          t = abs(brts[t] - brts[t - 1])
        ))
      }

      #Applying C operator (this is a trick to avoid precision issues)
      C[t] <- 1 / (sum(q_t[t, ])); q_t[t, ] <- q_t[t, ] * C[t]

      #Applying B operator
      if (t < max_t) {
        if (brts[t] != t_d) {
          q_t[t, ] <- q_t[t, ] * lambda
          k <- k + 1
        } else {
          q_t[t, ] <- q_t[t, ] #* (k + nvec)^-1 #NODIV
          k <- k - 1
        }

        #Applying D operator (this works exactly like C)
        D[t] <- 1 / (sum(q_t[t, ])); q_t[t, ] <- q_t[t, ] * D[t]

        #Updating running parameters
        t <- t + 1
      } else {
        break
      }
    }

    #Selecting the state I am interested in
    vm <- choose(k + missnumspec[clade], k) ^ -1
    p_m  <- vm * q_t[t, (missnumspec[clade] + 1)]

    #Removing C and D effects from the LL
    loglik <- log(p_m) - sum(log(C)) - sum(log(D))

    #Various checks
    loglik <- as.numeric(loglik)
    if (is.nan(loglik) | is.na(loglik)) {
      loglik <- -Inf
    }
    logliks[clade] <- loglik
  }

  pc <- sls::pc_1shift(
    pars_m = pars_m,
    pars_s = pars_s,
    brts_m = brts_m,
    brts_s = brts_s,
    cond = cond,
    n_max = n_max,
    n_0 = n_0
  )

  total_loglik <- sum(logliks) - log(pc)
  total_loglik <- as.numeric(total_loglik)
  if (is.nan(total_loglik) | is.na(total_loglik)) {
    total_loglik <- -Inf
  }
  return(total_loglik)
}

#' @title DDD-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood calling the routine from
#'  the DDD package
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_ddd <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  n_max = 1e2
) {
  pars_m <- pars[1:2]
  pars_s <- pars[3:4]
  brts_m <- brts[[1]]
  brts_s <- brts[[2]]
  pars1 <- c(pars_m[1], pars_m[2], Inf, pars_s[1], pars_s[2], Inf, brts_s[1])
  pars2 <- c(
    n_max,  #maximum number of species involved in the computation
    1,  #ddmodel: not actually used by this function
    0,  #conditioning
    min(abs(brts_m[abs(brts_m) > pars1[7]])), # tshift
    0,  #print things: not actually used by this function
    n_0 #stem or crown (soc)
  )
  brts_m1 <- brts_m
  brts_s1 <- brts_s[-1]
  missnumspec <- c(0, 0)

  testit::assert(pars2[4] %in% abs(brts_m1)) #tsplit in brts_m

  loglik <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = abs(brts_m1),
    brtsS = abs(brts_s1),
    missnumspec = missnumspec
  ); loglik

  pc <- sls::pc_1shift(
    pars_m = pars_m,
    pars_s = pars_s,
    brts_m = brts_m,
    brts_s = brts_s,
    cond = cond,
    n_max = n_max,
    n_0 = n_0
  ); pc

  loglik_c <- loglik - log(pc)
  return(loglik_c)
}

#' @title BISSE loglik with shift
#' @author Giovanni Laudanno
#' @description Provides BISSE loglik function in a presence of a shift.
#' @inheritParams default_params_doc
#' @return loglik
#' @export
loglik_bisse_shift <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  n_max = 1e2
) {
  pars_m <- pars[1:2]
  pars_s <- pars[3:4]
  brts_m <- brts[[1]]
  brts_s <- brts[[2]]
  loglik_s <- sls::loglik_bisse(
    pars = pars_s,
    brts = brts_s,
    n_0 = 1
  )
  loglik <- sls::loglik_bisse(
    pars = pars_m,
    brts = brts_m,
    n_0 = n_0,
    t_ds = brts_s[1],
    d_0s = exp(loglik_s)
  )

  pc <- sls::pc_1shift(
    pars_m = pars_m,
    pars_s = pars_s,
    brts_m = brts_m,
    brts_s = brts_s,
    cond = cond,
    n_max = n_max,
    n_0 = n_0
  ); pc

  loglik_c <- loglik - log(pc)
  return(loglik_c)
}

#' @title BISSE loglik shift
#' @author Giovanni Laudanno
#' @description Provides BISSE loglik shift function (alternative version).
#'  Yields the old (wrong) BISSE result for the Main Clade only.
#' @inheritParams default_params_doc
#' @return loglik
#' @export
loglik_bisse_shift2 <- function(
  pars,
  brts,
  n_0 = 2,
  t_0 = 0,
  t_d,
  log_scale = TRUE,
  lambdaterms = TRUE
) {
  testit::assert(all(t_d != brts))
  lambda <- pars[1]
  brts1 <- brts[brts > t_d]
  brts2 <- sort(c(t_d, brts[brts < t_d]), decreasing = TRUE)
  dd_1 <- sls::loglik_bisse2(
    pars,
    brts1,
    n_0 = n_0,
    t_0 = t_d,
    e_0 = sls::e_t(
      pars = pars,
      t_0 = t_0,
      t_f = t_d,
      e_0 = 0,
      d_0 = 1
    ),
    log_scale = FALSE,
    lambdaterms = FALSE
  )
  dd_2 <- sls::loglik_bisse2(
    pars,
    brts2,
    n_0 = (n_0 + length(brts1) - 1) - 1,
    t_0 = t_0,
    log_scale = FALSE,
    lambdaterms = FALSE
  )

  dd <- dd_1 * dd_2 * lambda ^ (length(brts[-1]) * lambdaterms)
  out <- (log_scale) * log(dd) + (1 - log_scale) * dd
  return(out)
}
