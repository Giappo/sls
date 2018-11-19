# MAIN COMPONENTS ----

#' @export
sls_sim.get_pars <- function(
  lambdas,
  mus,
  Ks = c(Inf, Inf)
) {
  n_clades <- length(lambdas)
  testit::assert(n_clades > 0)
  testit::assert(length(lambdas) == n_clades)
  testit::assert(length(mus) == n_clades)
  testit::assert(length(Ks) == n_clades)

  testit::assert(all(lambdas >= 0))
  testit::assert(all(mus >= 0))
  testit::assert(all(Ks > 0))

  pars <- vector("list", n_clades)
  for (clade in 1:n_clades) {
    pars[[clade]] <- c(
      lambda = lambdas[clade],
      mu = mus[clade],
      K = Ks[clade]
    )
    names(pars[[clade]]) <- c(
      "lambda",
      "mu",
      "K"
    )
  }
  pars
}

#' @export
sls_sim.initialize_LL_new_clade <- function(
  data,
  clade,
  pars,
  LS = sls:::sls_sim.get_standard_LS(),
  l_matrix_size = 1e4
) {
  if (clade == 0) {
    return(
      list(
        LL = vector("list", length(LS$clade_id)),
        pools = vector("list", length(LS$clade_id)),
        Nmax = vector("list", length(LS$clade_id))
      )
    )
  }

  if (l_matrix_size < 10) {
    l_matrix_size <- 10
  }

  LL <- data$LL
  n_0 <- sls_sim.get_n_0(LS = LS, clade = clade)
  t0 <- sls_sim.get_t0(LS = LS, clade = clade)
  motherclade <- sls_sim.get_motherclade(LS = LS, clade = clade)
  Nmax <- n_0

  cladeborn <- 1
  if (clade > 1) {
    Lmother <- LL[[motherclade]]
    motherspecies <- Lmother[Lmother[, 5] == clade, 3]
    cladeborn <- length(motherspecies) != 0
  }

  if (cladeborn) {
    t <- t0
    L <- matrix(0, nrow = l_matrix_size, 5)
    L[, 5] <- 0
    L[, 4] <- -1
    L[, 3] <- 0
    L[1, 1:4] <- c(t, 0,  1, -1)
    if (n_0 == 2) {
      L[2, 1:4] <- c(t, 1, -2, -1)
    }
    if (clade > 1) {
      L[1, 3] <- sign(motherspecies)
    }
    colnames(L) <- c(
      "birth_time",
      "parent",
      "id",
      "death_time",
      "shifted_to"
    )
    pool <- L[1:n_0, 3]
    data$LL[[clade]] <- L
    data$pools[[clade]] <- pool
    data$Nmax[[clade]] <- Nmax
  } else {
    data$LL[clade] <- list(NULL)
    data$pools[clade] <- list(NULL)
    data$Nmax[clade] <- list(NULL)
  }

  return(
    data
  )
}

#' @export
sls_sim.initialize_t_new_clade <- function(
  data,
  clade
) {
  LL <- data$LL

  if (!is.null(LL[[clade]])) {
    t <- unname(LL[[clade]][1, 1]); t
  } else {
    t <- 0
  }
  return(t)
}

#' @export
sls_sim.sample_deltas <- function(
  data,
  clade,
  pars
) {
  LL <- data$LL
  pools <- data$pools
  pool <- pools[[clade]]

  PARS   <- pars[[clade]]
  lambda <- PARS["lambda"]
  mu     <- PARS["mu"]

  N <- length(pool)
  total_rate <- N * (lambda + mu)
  testit::assert(total_rate >= 0)
  if (total_rate > 0) {
  delta_t <- (total_rate > 0) *
    rexp(1, rate = total_rate + (total_rate == 0)) +
    (total_rate == 0) * LL[[1]][1, 1]
  delta_n <- sample(c(-1, 1), size = 1, prob = c(mu, lambda))
  } else {
    delta_t <- 1000
    delta_n <- 0
  }
  return(list(
    delta_n = unname(delta_n),
    delta_t = unname(delta_t)
  ))
}

#' @export
sls_sim.decide_event <- function(
  data,
  clade,
  delta_n,
  delta_t,
  t,
  LS
) {
  LL <- data$LL
  L <- LL[[clade]]
  Nmax <- data$Nmax[[clade]]; L[1:(3 + 1), ]
  already_shifted <- any(L[1:Nmax, 5] > 0)
  tshifts <- sls:::sls_sim.get_shifts_info(LS = LS, clade = clade)
  if (nrow(tshifts) > 1) {
    stop("Check the function if you want to implement more than 1 shift!")
  }
  if (occur_shift <- (nrow(tshifts) > 0)) {
    occur_shift <- occur_shift &
      t > 0 &
      t > tshifts[clade] &
      (t - delta_t) < tshifts[clade] &
      already_shifted == 0
  }
  if (occur_shift) {
    return("shift")
  }

  if ( (t - delta_t) < 0 ) {
    return("end")
  }

  occur_speciation <- delta_n > 0 & (t - delta_t) > 0
  if (occur_speciation) {
    return("speciation")
  }

  occur_extinction <- delta_n < 0 & (t - delta_t) > 0
  if (occur_extinction) {
    return("extinction")
  }

  return("end")
}

#' @export
sls_sim.use_event <- function(
  data,
  clade,
  LS,
  event,
  t
) {

  shifts <- sls_sim.get_shifts_info(LS = LS, clade = clade)

  L <- data$LL[[clade]]
  pool <- data$pools[[clade]]; pool
  N <- length(pool)
  Nmax <- data$Nmax[[clade]]
  shifted <- 0

  if (event == "shift") {
    where <- shifts$where
    t <- shifts$when #time becomes the shift point
    if (N > 1) {
      shifted <- sample(pool, replace = FALSE, size = 1)
    } else {
      shifted <- pool
    }
    L[abs(shifted), 4] <- t #remove the shifted from the L table
    pool <- pool[pool != shifted]
    L[abs(L[, 3]) == abs(shifted), 5] <- where #register that the shift occurred
  }

  if (event == "speciation") {
    data <- sls_sim.adapt_l_matrix_size(
      data = data,
      clade = clade
    )
    L <- data$LL[[clade]]
    pool <- data$pools[[clade]]; pool
    N <- length(pool)

    if (N > 1) {
      parents <- sample(pool, replace = FALSE, size = 1)
    } else {
      parents <- pool
    }
    Nmax <- Nmax + 1

    new_line <- c(
      t,
      parents,
      abs(Nmax) * sign(parents),
      -1,
      0
    )
    dim(new_line) <- c(1, 5)
    L[Nmax, ] <- new_line
    pool <- c(pool, abs(Nmax) * sign(parents))
  }

  if (event == "extinction") {
    if (N > 1) {
      dead <- sample(pool, replace = FALSE, size = 1)
    } else {
      dead <- pool
    }
    L[abs(dead), 4] <- t
    pool <- pool[pool != dead]
  }

  if (event == "end" | length(pool) == 0) {
    t <- 0
    L2 <- L[1:Nmax, ]
    L <- L2
  }

  data2 <- data
  data2$LL[[clade]] <- L
  data2$pools[[clade]] <- pool
  data2$Nmax[[clade]] <- Nmax

  return(list(
    t = unname(t),
    data = data2
  ))
}

#' @export
sls_sim.check_survival <- function(
  L,
  final_time = 0
) {
  if (is.matrix(L)) {
    cond <- any(L[, 4] < final_time)
  } else {
    cond <- L[4] < final_time
  }
  return(cond)
}

#' @export
sls_sim.conditioning <- function(
  data,
  LS,
  cond
) {

  if (nrow(LS) > 2) {
    stop("Currently this only works for 1 shift!")
  }
  if (cond == 2) {
    stop("Cond 2 is not supported!")
  }

  LL <- data$LL; clade <- 1
  for (clade in LS$clade_id) {
    L <- LL[[clade]]
    if (!is.matrix(L) && !is.null(L)) {
      dim(L) <- c(1, 5)
    }
    shifts <- sls_sim.get_shifts_info(LS = LS, clade = clade)
    shifts_times <- shifts$when
    shifted_id <- L[(L[, 5] != 0), 3]

    tp <- 0; ts <- shifts_times[1]
    if (clade == 1) {
      # subclade does NOT start from here (M1)
      coords_left  <- sign(L[, 3]) == -sign(shifted_id)

      #subclade does start from here!!! (M2)
      coords_right <- sign(L[, 3]) ==  sign(shifted_id)

      L_left  <- L[coords_left, ]; dim(L_left)  <- c(sum(coords_left), 5)
      L_right <- L[coords_right, ]; dim(L_right) <- c(sum(coords_right), 5)
      if (nrow(L_right) > 0) {
        testit::assert(
          any(L_right[, 5] > 0)
        )
      }
      colnames(L_left) <- colnames(L_right) <- colnames(L)

      # lineages in M1 born before tshift
      coords_left_cs  <- (L_left[, 1]  > shifts_times[1])
      L_left_cs  <- L_left[coords_left_cs, ]
      dim(L_left_cs)  <- c(sum(coords_left_cs), 5)

      # lineages in M2 born before tshift
      coords_right_cs <- (L_right[, 1] > shifts_times[1])
      L_right_cs <- L_right[coords_right_cs, ]
      dim(L_right_cs) <- c(sum(coords_right_cs), 5)

      surv_left_cp  <- sls:::sls_sim.check_survival(
        L = L_left,
        final_time = tp
      ) #M1 survives from c to p
      surv_right_cp <- sls:::sls_sim.check_survival(
        L = L_right,
        final_time = tp
      ) #M2 survives from c to p
      surv_left_cs  <- sls:::sls_sim.check_survival(
        L = L_left_cs,
        final_time = ts
      ) #M1 survives from c to s
      surv_right_cs <- sls:::sls_sim.check_survival(
        L = L_right_cs,
        final_time = ts
      ) #M2 survives from c to s

      testit::assert(surv_left_cs  >= surv_left_cp)
      testit::assert(surv_right_cs >= surv_right_cp)
    } else {
      L <- LL[[clade]]
      surv_s <- sls:::sls_sim.check_survival(
        L = L,
        final_time = tp
      ) #S survives from s to p
    }
  }

  cond1 <- (surv_left_cp && surv_right_cs && surv_s) ||
    (surv_left_cp && surv_right_cp)
  cond2 <- 1
  cond3 <- (surv_left_cp && surv_right_cs && surv_s)
  cond4 <- (surv_left_cp && surv_right_cp && surv_s)

  #conditioning
  keep_the_sim <- (cond == 0) * 1 +
    (cond == 1) * cond1 +
    (cond == 2) * cond2 +
    (cond == 3) * cond3 +
    (cond == 4) * cond4
}

#' @export
sls_sim.get_brts <- function(
  data,
  LS
) {
  brts <- vector("list", length(data$LL))

  for (clade in seq_along(data$LL)) {
    done <- 0
    if (is.null(data$LL[[clade]]) && done == 0) {
      brts[clade] <- NULL
      done <- 1
    }
    L <- sls_sim.cut_l_matrix(
      unname(data$LL[[clade]])
    )
    if (!any(L[, 4] == -1) && done == 0) {
      brts[clade] <- NULL
      done <- 1
    }
    if (LS$n_0[clade] == 1 && done == 0) {
      if (sum(L[, 4] == -1) == 1 | nrow(L) == 1) {
        brts[[clade]] <- L[1, 1]
      } else {
        phylo <- DDD:::L2phylo(L, dropextinct = TRUE)
        brts[[clade]] <- ape::branching.times(phylo)
      }
    }
    if (LS$n_0[clade] == 2 && done == 0) {
      brts[[clade]] <- DDD:::L2brts(
        L,
        dropextinct = TRUE
      )
    }
  }

  brts
}

# UTILITIES ----

#' @export
sls_sim.read_LL <- function(
  LL
) {
  if (!is.list(LL)) {
    stop("LL is not a list!!!")
  }
  LL2 <- LL
  for (i in seq_along(LL)) {
    L <- LL[[i]]
    LL2[[i]] <- L[L[, 3] != 0, ]
  }
  LL2
}

#' @export
sls_sim.get_standard_LS <- function(
  crown_age = 10,
  n_0 = 2,
  shift_time = 4
) {
  testit::assert(crown_age > shift_time)
  LS <- as.data.frame(matrix(0, nrow = 2, ncol = 4))
  LS[, 1] <- c(0, shift_time)
  LS[, 2] <- c(0, 1)
  LS[, 3] <- c(1, 2)
  LS[, 4] <- c(n_0, 1)
  LS[1, 1] <- crown_age
  colnames(LS) <- c("birth_time", "motherclade", "clade_id", "n_0")

  t0s <- LS[, 1]
  motherclades <- LS[, 2]
  n_0s <- LS[, 4]
  testit::assert(
    all(motherclades < seq(from = 1, to = length(motherclades)))
  )
  testit::assert(all(n_0s >= 1))
  if (any(n_0s[-1] != 1)) {
    stop("Every subclade should start with 1 species!")
  }
  testit::assert(all(t0s > 0))

  LS
}

#' @export
sls_sim.get_shifts_info <- function(
  LS = sls:::sls_sim.get_standard_LS(),
  clade
) {
  when  <- LS[LS[, 2] == clade, 1]
  where <- LS[LS[, 2] == clade, 3]

  info <- data.frame(when = when, where = where)
  info
}

#' @export
sls_sim.get_t0 <- function(
  LS = sls:::sls_sim.get_standard_LS(),
  clade
) {
  t0 <- LS[LS[, 3] == clade, 1]
  t0
}

#' @export
sls_sim.get_n_0 <- function(
  LS = sls:::sls_sim.get_standard_LS(),
  clade
) {
  n_0 <- LS[LS[, 3] == clade, 4]
  n_0
}

#' @export
sls_sim.get_motherclade <- function(
  LS = sls:::sls_sim.get_standard_LS(),
  clade
) {
  motherclade <- LS[LS[, 3] == clade, 2]
  motherclade
}

#' @export
sls_sim.get_pool <- function(
  L
) {
  right_rows <- L[, 3] != 0 & L[, 4] == -1
  pool <- L[right_rows, 3]
  pool
}

#' @export
sls_sim.cut_l_matrix <- function(
  L
) {
  if (is.matrix(L)) {
    Nmax <- max(abs(L[, 3]))
    L2 <- L[1:Nmax, ]
    cols_L <- ncol(L)
  } else {
    Nmax <- max(abs(L[3]))
    L2 <- L
    cols_L <- length(L)
  }
  dim(L2) <- c(Nmax, cols_L)
  return(L2)
}

#' @export
sls_sim.adapt_l_matrix_size <- function(
  data,
  clade
) {
  Nmax <- data$Nmax[[clade]]
  L <- data$LL[[clade]]
  if (Nmax >= nrow(L) - 2) {
    append_L <- matrix(0, nrow = nrow(L), ncol = ncol(L))
    append_L[, 4] <- -1
    L2 <- rbind(L, append_L)
    data$LL[[clade]] <- L2
  }
  return(data)
}
