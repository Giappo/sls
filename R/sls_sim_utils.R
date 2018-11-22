# MAIN COMPONENTS ----

#' @title Get sls parameters
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return sls parameters
#' @export
sls_sim_get_pars <- function(
  lambdas,
  mus,
  ks = c(Inf, Inf)
) {
  n_clades <- length(lambdas)
  testit::assert(n_clades > 0)
  testit::assert(length(lambdas) == n_clades)
  testit::assert(length(mus) == n_clades)
  testit::assert(length(ks) == n_clades)

  testit::assert(all(lambdas >= 0))
  testit::assert(all(mus >= 0))
  testit::assert(all(ks > 0))

  pars <- vector("list", n_clades)
  for (clade in 1:n_clades) {
    pars[[clade]] <- c(
      lambda = lambdas[clade],
      mu = mus[clade],
      K = ks[clade]
    )
    names(pars[[clade]]) <- c(
      "lambda",
      "mu",
      "K"
    )
  }
  pars
}

#' @title Initialize data for a new clade
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return data for the clade
#' @export
sls_sim_initialize_data_new_clade <- function(
  data,
  clade,
  pars,
  l_2 = sls::sls_sim_get_standard_l_2(),
  l_matrix_size = 1e4
) {
  if (clade == 0) {
    return(
      list(
        l_1 = vector("list", length(l_2$clade_id)),
        pools = vector("list", length(l_2$clade_id)),
        n_max = vector("list", length(l_2$clade_id)),
        t = vector("list", length(l_2$clade_id))
      )
    )
  }

  if (l_matrix_size < 10) {
    l_matrix_size <- 10
  }

  l_1 <- data$l_1
  n_0 <- sls_sim_get_n_0(l_2 = l_2, clade = clade)
  t_0 <- sls_sim_get_t_0(l_2 = l_2, clade = clade)
  motherclade <- sls_sim_get_motherclade(l_2 = l_2, clade = clade)
  n_max <- n_0

  cladeborn <- 1
  if (clade > 1) {
    l_mother <- l_1[[motherclade]]
    motherspecies <- l_mother[l_mother[, 5] == clade, 3]
    cladeborn <- length(motherspecies) != 0
  }

  if (cladeborn) {
    t <- t_0
    l_0 <- matrix(0, nrow = l_matrix_size, 5)
    l_0[, 5] <- 0
    l_0[, 4] <- -1
    l_0[, 3] <- 0
    l_0[1, 1:4] <- c(t, 0,  1, -1)
    if (n_0 == 2) {
      l_0[2, 1:4] <- c(t, 1, -2, -1)
    }
    if (clade > 1) {
      l_0[1, 3] <- sign(motherspecies)
    }
    colnames(l_0) <- c(
      "birth_time",
      "parent",
      "id",
      "death_time",
      "shifted_to"
    )
    pool <- l_0[1:n_0, 3]
    data$l_1[[clade]] <- l_0
    data$pools[[clade]] <- pool
    data$n_max[[clade]] <- n_max
  } else {
    data$l_1[clade] <- list(NULL)
    data$pools[clade] <- list(NULL)
    data$n_max[clade] <- list(NULL)
  }

  data$t[[clade]] <- t_0
  return(
    data
  )
}

#' @title Sample deltas for Doob-Gillespie
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return delta_n and delta_t
#' @export
sls_sim_sample_deltas <- function(
  data,
  clade,
  pars
) {
  l_1 <- data$l_1
  pools <- data$pools
  pool <- pools[[clade]]

  PARS   <- pars[[clade]]
  lambda <- PARS["lambda"]
  mu     <- PARS["mu"]

  n <- length(pool)
  total_rate <- n * (lambda + mu)
  testit::assert(total_rate >= 0)
  if (total_rate > 0) {
  delta_t <- (total_rate > 0) *
    stats::rexp(1, rate = total_rate + (total_rate == 0)) +
    (total_rate == 0) * l_1[[1]][1, 1]
  delta_n <- sample(c(-1, 1), size = 1, prob = c(mu, lambda))
  } else {
    delta_t <- 1e5
    delta_n <- 0
  }
  return(list(
    delta_n = unname(delta_n),
    delta_t = unname(delta_t)
  ))
}

#' @title Determines the event
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return the event
#' @export
sls_sim_decide_event <- function(
  data,
  clade,
  deltas,
  l_2
) {
  delta_n <- deltas$delta_n
  delta_t <- deltas$delta_t
  t <- data$t[[clade]]
  l_1 <- data$l_1
  l_0 <- l_1[[clade]]
  already_shifted <- any(l_0[, 5] > 0)
  tshifts <- sls::sls_sim_get_shifts_info(l_2 = l_2, clade = clade)
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

#' @title Update data and time given the event
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return data
#' @export
sls_sim_use_event <- function(
  data,
  clade,
  l_2,
  event,
  deltas
) {

  shifts <- sls_sim_get_shifts_info(l_2 = l_2, clade = clade)
  t <- data$t[[clade]] - deltas$delta_t
  l_0 <- data$l_1[[clade]]
  pool <- data$pools[[clade]]; pool
  n <- length(pool)
  n_max <- data$n_max[[clade]]
  shifted <- 0

  if (event == "shift") {
    where <- shifts$where
    t <- shifts$when #time becomes the shift point
    if (n > 1) {
      shifted <- sample(pool, replace = FALSE, size = 1)
    } else {
      shifted <- pool
    }
    l_0[abs(shifted), 4] <- t #remove the shifted from the l_0 table
    pool <- pool[pool != shifted]
    l_0[abs(l_0[, 3]) == abs(shifted), 5] <- where #register if shift occurs
  }

  if (event == "speciation") {
    data <- sls_sim_adapt_l_matrix_size(
      data = data,
      clade = clade
    )
    l_0 <- data$l_1[[clade]]
    pool <- data$pools[[clade]]; pool
    n <- length(pool)

    if (n > 1) {
      parents <- sample(pool, replace = FALSE, size = 1)
    } else {
      parents <- pool
    }
    n_max <- n_max + 1

    new_line <- c(
      t,
      parents,
      abs(n_max) * sign(parents),
      -1,
      0
    )
    dim(new_line) <- c(1, 5)
    l_0[n_max, ] <- new_line
    pool <- c(pool, abs(n_max) * sign(parents))
  }

  if (event == "extinction") {
    if (n > 1) {
      dead <- sample(pool, replace = FALSE, size = 1)
    } else {
      dead <- pool
    }
    l_0[abs(dead), 4] <- t
    pool <- pool[pool != dead]
  }

  if (event == "end" | length(pool) == 0) {
    t <- 0
    l_02 <- sls_sim_cut_l_matrix(l_0)
    l_0 <- l_02
  }

  # store output
  t <- unname(t)
  data2 <- data
  if (is.null(l_0)) {
    data2$l_1[clade] <- list(l_0)
  } else {
    data2$l_1[[clade]] <- l_0
  }
  if (is.null(pool)) {
    data2$pools[clade] <- list(pool)
  } else {
    data2$pools[[clade]] <- pool
  }
  if (is.null(pool)) {
    data2$n_max[clade] <- list(n_max)
  } else {
    data2$n_max[[clade]] <- n_max
  }
  data2$t[[clade]] <- t

  return(
    data = data2
  )
}

#' @title Check if a clade survives untile final_time
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a boolean
#' @export
sls_sim_check_survival <- function(
  l_0,
  final_time = 0
) {
  if (is.matrix(l_0)) {
    cond <- any(l_0[, 4] < final_time)
  } else {
    cond <- l_0[4] < final_time
  }
  return(cond)
}

#' @title Check if the conditioning is fulfilled
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a boolean
#' @export
sls_sim_conditioning <- function(
  data,
  l_2,
  cond
) {

  if (nrow(l_2) > 2) {
    stop("Currently this only works for 1 shift!")
  }
  if (cond == 2) {
    stop("Cond 2 is not supported!")
  }

  l_1 <- data$l_1; clade <- 1
  for (clade in l_2$clade_id) {
    l_0 <- l_1[[clade]]
    if (!is.matrix(l_0) && !is.null(l_0)) {
      dim(l_0) <- c(1, 5)
    }
    shifts <- sls_sim_get_shifts_info(l_2 = l_2, clade = clade)
    shifts_times <- shifts$when
    shifted_id <- l_0[(l_0[, 5] != 0), 3]

    t_p <- 0; ts <- shifts_times[1]
    if (clade == 1) {
      # subclade does NOT start from here (M1)
      coords_left  <- sign(l_0[, 3]) == -sign(shifted_id)

      #subclade does start from here!!! (M2)
      coords_right <- sign(l_0[, 3]) ==  sign(shifted_id)

      l_0_left  <- l_0[coords_left, ]
      dim(l_0_left)  <- c(sum(coords_left), 5)
      l_0_right <- l_0[coords_right, ]
      dim(l_0_right) <- c(sum(coords_right), 5)
      if (nrow(l_0_right) > 0) {
        testit::assert(
          any(l_0_right[, 5] > 0)
        )
      }
      colnames(l_0_left) <- colnames(l_0_right) <- colnames(l_0)

      # lineages in M1 born before tshift
      coords_left_cs  <- (l_0_left[, 1]  > shifts_times[1])
      l_0_left_cs  <- l_0_left[coords_left_cs, ]
      dim(l_0_left_cs)  <- c(sum(coords_left_cs), 5)

      # lineages in M2 born before tshift
      coords_right_cs <- (l_0_right[, 1] > shifts_times[1])
      l_0_right_cs <- l_0_right[coords_right_cs, ]
      dim(l_0_right_cs) <- c(sum(coords_right_cs), 5)

      surv_left_cp  <- sls::sls_sim_check_survival(
        l_0 = l_0_left,
        final_time = t_p
      ) #M1 survives from c to p
      surv_right_cp <- sls::sls_sim_check_survival(
        l_0 = l_0_right,
        final_time = t_p
      ) #M2 survives from c to p
      surv_left_cs  <- sls::sls_sim_check_survival(
        l_0 = l_0_left_cs,
        final_time = ts
      ) #M1 survives from c to s
      surv_right_cs <- sls::sls_sim_check_survival(
        l_0 = l_0_right_cs,
        final_time = ts
      ) #M2 survives from c to s

      testit::assert(surv_left_cs  >= surv_left_cp)
      testit::assert(surv_right_cs >= surv_right_cp)
    } else {
      l_0 <- l_1[[clade]]
      surv_s <- sls::sls_sim_check_survival(
        l_0 = l_0,
        final_time = t_p
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

#' @title Extract the branching times from data
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return the branching times
#' @export
sls_sim_get_brts <- function(
  data,
  l_2
) {
  if (!all(data$t == 0)) {
    stop("times in all clades need to reach 0 before you call this function")
  }
  brts <- vector("list", length(data$l_1))

  for (clade in seq_along(data$l_1)) {
    done <- 0
    if (is.null(data$l_1[[clade]]) && done == 0) {
      brts[clade] <- NULL
      done <- 1
    }
    if (done == 0) {
      l_0 <- sls_sim_cut_l_matrix(
        unname(data$l_1[[clade]])
      )
    }
    if (!any(l_0[, 4] == -1) && done == 0) {
      brts[clade] <- NULL
      done <- 1
    }
    if (l_2$n_0[clade] == 1 && done == 0) {
      if (sum(l_0[, 4] == -1) == 1 | nrow(l_0) == 1) {
        brts[[clade]] <- l_0[1, 1]
      } else {
        phylo <- DDD::L2phylo(l_0, dropextinct = TRUE)
        brts[[clade]] <- c(
          l_0[1, 1],
          ape::branching.times(phylo)
        )
      }
    }
    if (l_2$n_0[clade] == 2 && done == 0) {
      brts[[clade]] <- DDD::L2brts(
        l_0,
        dropextinct = TRUE
      )
    }
  }
  brts <- unname(brts)
  brts
}

# UTILITIES ----

#' @title Reads l_1
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return l_1
#' @export
sls_sim_read_l_1 <- function(
  l_1
) {
  if (!is.list(l_1)) {
    stop("l_1 is not a list!!!")
  }
  l_12 <- l_1
  for (i in seq_along(l_1)) {
    l_0 <- l_1[[i]]
    l_12[[i]] <- l_0[l_0[, 3] != 0, ]
  }
  l_12
}

#' @title Creates the standard l_2
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return l_2
#' @export
sls_sim_get_standard_l_2 <- function(
  crown_age = 10,
  n_0 = 2,
  shift_time = 4
) {
  testit::assert(crown_age > shift_time)
  l_2 <- as.data.frame(matrix(0, nrow = 2, ncol = 4))
  l_2[, 1] <- c(0, shift_time)
  l_2[, 2] <- c(0, 1)
  l_2[, 3] <- c(1, 2)
  l_2[, 4] <- c(n_0, 1)
  l_2[1, 1] <- crown_age
  colnames(l_2) <- c("birth_time", "motherclade", "clade_id", "n_0")

  t_0s <- l_2[, 1]
  motherclades <- l_2[, 2]
  n_0s <- l_2[, 4]
  testit::assert(
    all(motherclades < seq(from = 1, to = length(motherclades)))
  )
  testit::assert(all(n_0s >= 1))
  if (any(n_0s[-1] != 1)) {
    stop("Every subclade should start with 1 species!")
  }
  testit::assert(all(t_0s > 0))

  l_2
}

#' @title Get shifts' information
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return when and where the shifts occurs
#' @export
sls_sim_get_shifts_info <- function(
  l_2 = sls::sls_sim_get_standard_l_2(),
  clade
) {
  when  <- l_2[l_2[, 2] == clade, 1]
  where <- l_2[l_2[, 2] == clade, 3]

  info <- data.frame(when = when, where = where)
  info
}

#' @title Initialize t for a new clade
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return t_0 for the clade
#' @export
sls_sim_get_t_0 <- function(
  l_2 = sls::sls_sim_get_standard_l_2(),
  clade
) {
  t_0 <- l_2[l_2[, 3] == clade, 1]
  t_0
}

#' @title Initialize n for a new clade
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return n0 for the clade
#' @export
sls_sim_get_n_0 <- function(
  l_2 = sls::sls_sim_get_standard_l_2(),
  clade
) {
  n_0 <- l_2[l_2[, 3] == clade, 4]
  n_0
}

#' @title Get the motherclade of a given clade
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return the motherclade
#' @export
sls_sim_get_motherclade <- function(
  l_2 = sls::sls_sim_get_standard_l_2(),
  clade
) {
  motherclade <- l_2[l_2[, 3] == clade, 2]
  motherclade
}

#' @title Get the pool from the l_0 table
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return the pool of species
#' @export
sls_sim_get_pool <- function(
  l_0
) {
  right_rows <- l_0[, 3] != 0 & l_0[, 4] == -1
  pool <- l_0[right_rows, 3]
  pool
}

#' @title Cuts the l_0 table only to the useful rows
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a clean l_0 table
#' @export
sls_sim_cut_l_matrix <- function(
  l_0
) {
  if (is.null(l_0)) {
    return(NULL)
  }
  if (is.matrix(l_0)) {
    n_max <- max(abs(l_0[, 3]))
    l_02 <- l_0[1:n_max, ]
    cols_l_0 <- ncol(l_0)
  } else {
    n_max <- max(abs(l_0[3]))
    l_02 <- l_0
    cols_l_0 <- length(l_0)
  }
  dim(l_02) <- c(n_max, cols_l_0)
  return(l_02)
}

#' @title Increase the size of the l_0 table if needed
#' @description sls_sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a bigger l_0 table
#' @export
sls_sim_adapt_l_matrix_size <- function(
  data,
  clade
) {
  n_max <- data$n_max[[clade]]
  l_0 <- data$l_1[[clade]]
  if (n_max >= nrow(l_0) - 2) {
    append_l_0 <- matrix(0, nrow = nrow(l_0), ncol = ncol(l_0))
    append_l_0[, 4] <- -1
    l_02 <- rbind(l_0, append_l_0)
    data$l_1[[clade]] <- l_02
  }
  return(data)
}
