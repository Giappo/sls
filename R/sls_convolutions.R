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
#' @description Convolves all the processes before the shift and
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
#' @description Convolves all the processes before the shift and
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
  dft_p_t_n <- matrix(
    unlist(lapply(p_t_n, FUN = sls::dft)),
    nrow = n_t,
    byrow = TRUE
  )
  rownames(dft_p_t_n) <- paste0("t", 1:n_t)
  Re(sum(sls::idft(apply(dft_p_t_n, MARGIN = 2, "prod")))) #awesome!
}


#' @title Combine pn
#' @author Giovanni Laudanno
#' @description Uses the function DDD::conv to convolve
#'  all the processes before the shift and imposes the death before the present
#'  of all species that are not visible in the phylogeny
#' @inheritParams default_params_doc
#' @return Convolution of the probabilities for all the processes
#' @export
combine_ddd <- function(
  lambda,
  mu,
  times,
  tbar,
  n_max = 1e2
) {
  nvec <- 0:n_max
  n_t <- length(times)
  p_t_n <- vector("list", n_t)
  for (t in 1:n_t) {
    p_t_n[[t]] <- sls::pn_bar(
      n = nvec,
      t = times[t],
      lambda = lambda,
      mu = mu,
      tbar = tbar
    )
  }
  if (n_t == 1) {
    convol <- p_t_n[[1]]
  }
  if (n_t == 2) {
    convol <- DDD::conv(p_t_n[[1]], p_t_n[[2]])
  }
  if (n_t > 2) {
    x <- p_t_n
    while (length(x) > 1) {
      n1 <- length(x)
      n2 <- ceiling(n1 / 2)
      convol <- list()
      for (i in 1:n2) {
        j <- 2 * i - 1
        if (j > n1 - 1) {
          convol[[i]] <- x[[j]]
        } else {
          convol[[i]] <- DDD::conv(x[[j]], x[[j + 1]])
        }
      }
      x <- convol
    }
    convol <- unlist(x)
  }
  nvec <- c(Inf, 1:(length(convol) - 1))
  out <- sum(
    (nvec ^ -1) * convol
  )
  out
}

#' @title Combine pn in the old wrong way
#' @author Giovanni Laudanno
#' @description Uses the function DDD::conv to convolve all the processes before
#'  the shift and imposes the death before the present of all species that
#'  are not visible in the phylogeny.
#'  It doesn't divide by N. Used to check on old models.
#' @inheritParams default_params_doc
#' @return Convolution of the probabilities for all the processes
#' @export
combine_ddd_nodiv <- function(
  lambda,
  mu,
  times,
  tbar,
  n_max = 1e2
) {
  nvec <- 0:n_max
  n_t <- length(times)
  p_t_n <- vector("list", n_t)
  for (t in 1:n_t) {
    p_t_n[[t]] <- sls::pn_bar(
      n = nvec,
      t = times[t],
      lambda = lambda,
      mu = mu,
      tbar = tbar
    )
  }
  if (n_t == 1) {
    convol <- p_t_n[[1]]
  }
  if (n_t == 2) {
    convol <- DDD::conv(p_t_n[[1]], p_t_n[[2]])
  }
  if (n_t > 2) {
    x <- p_t_n
    while (length(x) > 1) {
      n1 <- length(x)
      n2 <- ceiling(n1 / 2)
      convol <- list()
      for (i in 1:n2) {
        j <- 2 * i - 1
        if (j > n1 - 1) {
          convol[[i]] <- x[[j]]
        } else {
          convol[[i]] <- DDD::conv(x[[j]], x[[j + 1]])
        }
      }
      x <- convol
    }
    convol <- unlist(x)
  }
  nvec <- c(Inf, 1:(length(convol) - 1))
  out <- sum(convol)
  out
}

#' @title Calculates the preshift loglik using convolutions
#' @author Giovanni Laudanno
#' @description Calculates the preshift loglik using convolutions
#' @inheritParams default_params_doc
#' @return Loglik preshift for main clade
#' @export
loglik_preshift <- function(
  lambdas,
  mus,
  brts,
  n_0 = 2,
  n_max = 1e2
) {
  brts_m <- brts[[1]]
  brts_s <- brts[[2]]
  brts_m1 <- sort(brts_m, decreasing = TRUE)
  brts_s1 <- sort(brts_s, decreasing = TRUE)
  t_d <- brts_s1[1]
  if (n_0 == 2) {
    brts_m2 <- c(brts_m1[1], brts_m1)
  } else {
    brts_m2 <- brts_m1
  }

  if (is.list(brts)) {
    n_min <- 2 * (0 + (n_0 - 1) + length(unlist(brts[[1]])))
  } else {
    n_min <- 2 * (0 + (n_0 - 1) + length(brts))
  }
  if (n_max < n_min) {
    n_max <- n_min
  }

  delta_s_i <- brts_m2[brts_m2 > t_d] - t_d
  k_shift <- length(delta_s_i)
  delta_p_s <- t_d

  one_minus_p_p_s <-
    one_minus_pt(lambda = lambdas[1], mu = mus[1], t = delta_p_s)
  p_s_i <-
    p_t(lambda = lambdas[1], mu = mus[1], t = delta_s_i)
  one_minus_u_s_i <-
    one_minus_ut(lambda = lambdas[1], mu = mus[1], t = delta_s_i)
  u_s_i <-
    ut(lambda = lambdas[1], mu = mus[1], t = delta_s_i)
  names(one_minus_u_s_i) <- names(p_s_i) <- names(u_s_i) <-
    paste0("t=", 1:k_shift)

  const_term <- sum(log(p_s_i) + log(one_minus_u_s_i))

  nvec <- 1:n_max
  p_t_n <- vector("list", k_shift)
  for (t in 1:k_shift) {
    p_t_n[[t]] <- nvec * (u_s_i[t] ^ (nvec - 1))
  }
  convolution2 <- p_t_n[[1]]
  fixed_term <- 0
  for (i in 2:k_shift) {
    convolution2 <- DDD::conv(p_t_n[[i]], convolution2[1:n_max])
    norm_term <- sum(convolution2)
    convolution2 <- convolution2 / norm_term
    fixed_term <- fixed_term + log(norm_term)
  }
  convolution3 <- c(rep(0, k_shift - 1), convolution2)
  convolution4 <- convolution3[1:n_max]
  norm_term <- sum(convolution4)
  convolution <- convolution4 / norm_term
  fixed_term <- fixed_term + log(norm_term)

  vec_one_over_n <- (nvec ^ -1)[k_shift:n_max]
  vec_one_minus_p <- (one_minus_p_p_s ^ (nvec - k_shift))[k_shift:n_max]
  vec_convolution <- convolution[k_shift:n_max]

  find_exp <- function(n_max) {
    300 * (1 - exp(-n_max / 1000))
  }

  safe_term_1 <- 2 ^ -find_exp(n_max)
  vec_one_over_n <- vec_one_over_n / safe_term_1
  fixed_term <- fixed_term + log(safe_term_1)

  vec_one_minus_p <- vec_one_minus_p / safe_term_1
  fixed_term <- fixed_term + log(safe_term_1)

  sum_term <- sum(
    (vec_one_over_n * vec_one_minus_p * vec_convolution)
  ) # awesome!

  loglik <- log(k_shift) + const_term + log(sum_term) + fixed_term
  loglik
}

#' @title Calculates the preshift loglik using convolutions
#' @author Giovanni Laudanno
#' @description Calculates the preshift loglik using convolutions
#' @inheritParams default_params_doc
#' @return Loglik preshift for main clade
#' @export
loglik_preshift_nodiv <- function(
  lambdas,
  mus,
  brts,
  n_0 = 2,
  n_max = 1e2
) {
  brts_m <- brts[[1]]
  brts_s <- brts[[2]]
  brts_m1 <- sort(brts_m, decreasing = TRUE)
  brts_s1 <- sort(brts_s, decreasing = TRUE)
  t_d <- brts_s1[1]
  if (n_0 == 2) {
    brts_m2 <- c(brts_m1[1], brts_m1)
  } else {
    brts_m2 <- brts_m1
  }

  if (is.list(brts)) {
    n_min <- 2 * (0 + (n_0 - 1) + length(unlist(brts[[1]])))
  } else {
    n_min <- 2 * (0 + (n_0 - 1) + length(brts))
  }
  if (n_max < n_min) {
    n_max <- n_min
  }

  delta_s_i <- brts_m2[brts_m2 > t_d] - t_d
  k_shift <- length(delta_s_i)
  delta_p_s <- t_d

  one_minus_p_p_s <-
    one_minus_pt(lambda = lambdas[1], mu = mus[1], t = delta_p_s)
  p_s_i <-
    p_t(lambda = lambdas[1], mu = mus[1], t = delta_s_i)
  one_minus_u_s_i <-
    one_minus_ut(lambda = lambdas[1], mu = mus[1], t = delta_s_i)
  u_s_i <-
    ut(lambda = lambdas[1], mu = mus[1], t = delta_s_i)
  names(one_minus_u_s_i) <- names(p_s_i) <- names(u_s_i) <-
    paste0("t=", 1:k_shift)

  const_term <- sum(log(p_s_i) + log(one_minus_u_s_i))

  nvec <- 1:n_max
  p_t_n <- vector("list", k_shift)
  for (t in 1:k_shift) {
    p_t_n[[t]] <- nvec * (u_s_i[t] ^ (nvec - 1))
  }
  convolution2 <- p_t_n[[1]]
  fixed_term <- 0
  for (i in 2:k_shift) {
    convolution2 <- DDD::conv(p_t_n[[i]], convolution2[1:n_max])
    norm_term <- sum(convolution2)
    convolution2 <- convolution2 / norm_term
    fixed_term <- fixed_term + log(norm_term)
  }
  convolution3 <- c(rep(0, k_shift - 1), convolution2)
  convolution4 <- convolution3[1:n_max]
  norm_term <- sum(convolution4)
  convolution <- convolution4 / norm_term
  fixed_term <- fixed_term + log(norm_term)

  vec_one_minus_p <- (one_minus_p_p_s ^ (nvec - k_shift))[k_shift:n_max]
  vec_convolution <- convolution[k_shift:n_max]

  find_exp <- function(n_max) {
    300 * (1 - exp(-n_max / 1000))
  }

  safe_term_1 <- 2 ^ -find_exp(n_max)

  vec_one_minus_p <- vec_one_minus_p / safe_term_1
  fixed_term <- fixed_term + log(safe_term_1)

  sum_term <- sum(
    (vec_one_minus_p * vec_convolution)
  ) # awesome!

  loglik <- const_term + log(sum_term) + fixed_term
  loglik
}
