#' @title P-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convolving Nee's functions
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_sls_p <- function(
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
  if (any(c(pars_m, pars_s) < 0) | any(c(pars_m, pars_s) > 90)) {
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

  sls::sls_check_input(
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

  if (n_0 == 2) {
    brts_m2 <- c(brts_m1[1], brts_m1)
  } else {
    brts_m2 <- brts_m1
  }
  ts_m_pre_shift  <- brts_m2[brts_m2 > t_d] - t_d; ts_m_pre_shift
  ts_m_post_shift <- brts_m2[brts_m2 < t_d] ; ts_m_post_shift
  k_shift <- length(ts_m_pre_shift)
  if (length(ts_m_post_shift) == 0) {
    ts_m_post_shift <- 0
  }
  if (length(ts_m_pre_shift) == 0) {
    cat("There are no branching times before the shift"); return(-Inf)
  }

  log_lik_m_pre_shift <- loglik_preshift(
    lambdas = lambdas,
    mus = mus,
    brts = brts,
    n_0 = n_0,
    n_max = n_max
  )
  log_lik_m_post_shift <- sum(
    log(
      sls::pn(
        n = 1,
        lambda = lambdas[1],
        mu = mus[1],
        t = ts_m_post_shift
      )
    )
  ) +
    (k_shift - 1) *
    log(
      sls::pn(
        n = 1,
        t = t_d,
        lambda = lambdas[1],
        mu = mus[1]
      )
    )
  log_lik_s_post_shift <- sum(
    log(sls::pn(n = 1, lambda = lambdas[2], mu = mus[2], t = brts_s1))
  ); log_lik_s_post_shift
  loglik_m0 <- log_lik_m_pre_shift + log_lik_m_post_shift
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

#' @title P-likelihood
#' @author Giovanni Laudanno
#' @description Calculates the likelihood convolving Nee's functions
#'  It works in a less efficient way.
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_sls_p2 <- function(
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
  if (any(c(pars_m, pars_s) < 0) | any(c(pars_m, pars_s) > 90)) {
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
  events_matrix <- (mat <- cbind(brts_matrix, t_d_matrix))[, order(-mat[1, ])]
  # this is k vector after the events
  kvec_m_after <- (n_0 - 1) + cumsum(events_matrix[2, ])
  # this is k vector before the events
  kvec_m_before <- c(n_0 - 1, kvec_m_after[-length(kvec_m_after)])
  k_shift <- kvec_m_before[events_matrix[2, ] == -1]

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

  lik_m_pre_shift  <- k_shift *
    sls::combine_pns(
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
  ) +
    (length(ts_m_pre_shift) - 1) *
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
  ); log_lik_s_post_shift
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
#' @description Calculates the likelihood integrating the Q-equation
#' @inheritParams default_params_doc
#' @return The likelihood
#' @export
loglik_sls_q <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  n_max = 1e2,
  missnumspec = 0
) {
  pars_m <- pars[1:2]
  pars_s <- pars[3:4]
  brts_m <- brts[[1]]
  brts_s <- brts[[2]]
  n_shifts <- length(brts) - 1
  k_m <- n_0 + length(brts_m) - n_shifts
  k_s <- length(brts_s)
  if (any(c(pars_m, pars_s) < 0) | any(c(pars_m, pars_s) > 90)) {
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

  sls::sls_check_input(
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
  nvec <- 0:n_max

  missnumspec_vec <- 0:missnumspec
  logliks_vec <- rep(0, length(missnumspec_vec))
  for (m_m in missnumspec_vec) {
    m_s <- missnumspec - m_m
    missnumspecs <- c(m_m, m_s)
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
      sumvec2 <- sumvec1 <- rep(1, max_t)

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
          # q_t[t, ] <- abs(expoRkit::expv( # nolint
          #   v = q_t[(t - 1), ], # nolint
          #   x = transition_matrix, # nolint
          #   t = abs(brts[t] - brts[t - 1]) # nolint
          # )) # nolint
          transition_matrix <- as.matrix(transition_matrix)
          q_t[t, ] <- abs(deSolve::ode(
            y = q_t[(t - 1), ],
            parms = transition_matrix,
            times = c(0, abs(brts[t] - brts[t - 1])),
            func = sls::sls_loglik_rhs
          )[2, -1])
        }

        #Applying sumvec1 operator (this is a trick to avoid precision issues)
        sumvec1[t] <- 1 / (sum(q_t[t, ])); q_t[t, ] <- q_t[t, ] * sumvec1[t]

        #Applying B operator
        if (t < max_t) {
          if (brts[t] != t_d) {
            q_t[t, ] <- q_t[t, ] * lambda
            k <- k + 1
          } else {
            q_t[t, ] <- q_t[t, ] * k * (k + nvec) ^ -1
            k <- k - 1
          }

          #Applying sumvec2 operator (this works exactly like sumvec1)
          sumvec2[t] <- 1 / (sum(q_t[t, ])); q_t[t, ] <- q_t[t, ] * sumvec2[t]

          #Updating running parameters
          t <- t + 1
        } else {
          break
        }
      }

      #Selecting the state I am interested in
      vm <- choose(k + missnumspecs[clade], k) ^ -1
      p_m  <- vm * q_t[t, (missnumspecs[clade] + 1)]

      #Removing sumvec1 and sumvec2 effects from the LL
      loglik <- log(p_m) - sum(log(sumvec1)) - sum(log(sumvec2))

      #Various checks
      loglik <- as.numeric(loglik)
      if (is.nan(loglik) | is.na(loglik)) {
        loglik <- -Inf
      }
      logliks[clade] <- loglik
    }
    logliks_vec[m_m + 1] <- sum(logliks) + lchoose(missnumspec, m_m)
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

  # total_loglik <- sum(logliks) - log(pc)
  total_loglik <- log(sum(exp(logliks_vec))) - log(pc)
  total_loglik <- as.numeric(total_loglik)
  if (is.nan(total_loglik) | is.na(total_loglik)) {
    total_loglik <- -Inf
  }
  return(total_loglik)
}
