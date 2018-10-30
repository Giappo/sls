# MAIN SCRIPT ----
#' @export
sls_sim2 <- function(
  lambdas,
  mus,
  Ks = c(Inf, Inf),
  cond = 3,
  LS = sls:::sls_sim.get_standard_LS()
)
{
  # rm(list = ls()); lambdas <- c(0.2, 0.4); mus <- c(0.1, 0.05); Ks <- c(Inf, Inf); cond <- 3; LS <- sls_sim.get_standard_LS()
  pars <- sls_sim.get_pars(
    lambdas = lambdas,
    mus = mus,
    Ks = Ks
  )
  good_sim <- 0
  while (!good_sim)
  {
    LL <- sls_sim.initialize_LL_new_clade(clade = 0); clade <- 1
    for (clade in LS$clade_id)
    {
      LL <- sls_sim.initialize_LL_new_clade(
        LL = LL,
        clade = clade,
        pars = pars,
        LS = LS
      ); #sls_sim.read_LL(LL = LL)
      t <- sls_sim.initialize_t_new_clade(
        LL = LL,
        clade = clade
      )
      while (t > 0)
      {
        deltas <- sls_sim.sample_deltas(
          LL = LL,
          clade = clade,
          pars = pars
        ); deltaN <- deltas$deltaN; deltaT <- deltas$deltaT;
        event <- sls_sim.decide_event(
          LL = LL,
          clade = clade,
          LS = LS,
          deltaN = deltaN,
          deltaT = deltaT,
          t = t
        ); event
        output <- sls_sim.use_event(
          LL = LL,
          clade = clade,
          LS = LS,
          event = event,
          t = t - deltas$deltaT
        ); head(output$L); output$t
        t <- output$t
        LL[[clade]] <- output$L
      }; output$L
    }; LL
    good_sim <- sls_sim.conditioning(
      LL = LL,
      LS = LS,
      cond = cond
    ); good_sim
  }
  return(LL)
}

# MAIN COMPONENTS ----

#' @export
sls_sim.get_pars <- function(
  lambdas,
  mus,
  Ks = c(Inf, Inf)
)
{
  Nclades <- length(lambdas)
  testit::assert(Nclades > 0)
  testit::assert(length(lambdas) == Nclades)
  testit::assert(length(mus) == Nclades)
  testit::assert(length(Ks) == Nclades)

  testit::assert(all(lambdas >= 0))
  testit::assert(all(mus >= 0))
  testit::assert(all(Ks > 0))

  pars <- vector("list", Nclades)
  for (clade in 1:Nclades)
  {
    pars[[clade]] <- c(
      lambda = lambdas[clade],
      mu = mus[clade],
      K = Ks[clade]
    )
    names(pars[[clade]]) <- c("lambda",
                              "mu",
                              "K"
    )
  }
  pars
}

#' @export
sls_sim.initialize_LL_new_clade <- function(
  LL,
  clade,
  pars,
  LS = sls:::sls_sim.get_standard_LS()
)
{
  if (clade == 0) {return(LL = vector("list", length(LS$clade_id)))}

  PARS   <- pars[[clade]]
  lambda <- PARS['lambda']
  mu     <- PARS['mu']
  K      <- PARS['K']
  N0 <- sls_sim.get_N0(LS = LS, clade = clade)
  t0 <- sls_sim.get_t0(LS = LS, clade = clade)
  motherclade <- sls_sim.get_motherclade(LS = LS, clade = clade)

  cladeborn <- 1
  if (clade > 1)
  {
    Lmother <- LL[[motherclade]]
    motherspecies <- Lmother[Lmother[, 5] == clade, 3]
    cladeborn <- length(motherspecies) != 0
  }

  if (cladeborn) {
    total_count <- N0
    pool <- 1:N0
    N <- N0
    t <- t0
    L <- matrix(0, nrow = 1e6, 5)
    L[, 5] <- 0
    L[, 4] <- -1
    L[, 3] <- 0
    L[1, 1:4] <- c(t, 0,  1, -1)
    if (N0 == 2){L[2, 1:4] <- c(t, 1, -2, -1)}; head(L)
    if (clade > 1) {
      L[1, 3] <- sign(motherspecies)
    }
    colnames(L) <- c("birth_time",
                     "parent",
                     "id",
                     "death_time",
                     "shifted_to")
    LL[[clade]] <- L
  } else {
    LL[clade] <- list(NULL)
  }

  return(LL)
}

#' @export
sls_sim.initialize_t_new_clade <- function(
  LL,
  clade
)
{
  if (!is.null(LL[[clade]]))
  {
    t <- unname(LL[[clade]][1, 1]); t
  }
  else
  {
    t <- 0
  }
  return(t)
}

#' @export
sls_sim.sample_deltas <- function(
  LL,
  clade,
  pars
)
{
  L <- LL[[clade]]

  pool   <- sls_sim.get_pool(L)
  PARS   <- pars[[clade]]
  lambda <- PARS['lambda']
  mu     <- PARS['mu']

  N <- length(pool)
  total_rate <- N * (lambda + mu)
  deltaT <- (total_rate > 0) * rexp(1, rate = total_rate + ((total_rate == 0))) +
    (total_rate == 0) * LL[[1]][1, 1]
  deltaN <- sample(c(-1, 1), size = 1, prob = c(mu, lambda))
  return(list(
    deltaN = unname(deltaN),
    deltaT = unname(deltaT)
  ))
}

#' @export
sls_sim.decide_event <- function(
  LL,
  clade,
  deltaN,
  deltaT,
  t,
  LS
)
{
  L <- LL[[clade]]
  already_shifted <- sum(L[, 5]) > 0
  tshifts <- sls:::sls_sim.get_shifts_info(LS = LS, clade = clade)
  if (nrow(tshifts) > 1) {stop('Check the function if you want to implement more than 1 shift!')}

  if (occur_shift <- (nrow(tshifts) > 0))
  {
    occur_shift <- occur_shift &
      t > 0 &
      t > tshifts[clade] &
      (t - deltaT) < tshifts[clade] &
      already_shifted == 0
  }
  if (occur_shift) {return('shift')}

  if ((t - deltaT) < 0) {return('end')}

  occur_speciation <- deltaN > 0 & (t - deltaT) > 0
  if (occur_speciation) {return('speciation')}


  occur_extinction <- deltaN < 0 & (t - deltaT) > 0
  if (occur_extinction) {return('extinction')}

  return('end')
}

#' @export
sls_sim.use_event <- function(
  LL,
  clade,
  LS,
  event,
  t
)
{

  shifts  <- sls_sim.get_shifts_info(LS = LS, clade = clade)

  L <- LL[[clade]]
  pool <- sls_sim.get_pool(L)
  N <- length(pool)
  Nmax <- max(abs(L[,3]))
  shifted <- 0

  if (event == 'shift')
  {
    where <- shifts$where
    t <- shifts$when #time becomes the shift point
    if (N > 1) {shifted <- sample(pool, replace = FALSE, size = 1)}else{shifted <- pool}
    L[abs(shifted), 4] <- t #remove the shifted from the L table
    pool <- pool[pool != shifted]
    L[L[, 3] == shifted, 5] <- where #register that the shift occurred
  }

  if (event == 'speciation')
  {
    if (N > 1) {
      parents <- sample(pool, replace = FALSE, size = 1)}else{parents <- pool}
    new_line <- c(
      t,
      parents,
      abs(Nmax + 1) * sign(parents),
      -1,
      0
    )
    L[Nmax + 1, ] <- new_line
  }

  if (event == 'extinction')
  {
    if (N > 1) {dead <- sample(pool, replace = FALSE, size = 1)}else{dead <- pool}
    L[abs(dead), 4] <- t
  }

  if (event == 'end' | length(sls_sim.get_pool(L)) == 0)
  {
    t <- 0
    L2 <- sls_sim.cut_L(L)
    L <- L2
  }

  return(list(
    t = unname(t),
    L = L
  ))
}

#' @export
sls_sim.check_survival <- function(
  L,
  final_time = 0
)
{
  if (is.matrix(L))
  {
    cond <- any(L[,4] < final_time)
  }else
  {
    cond <- L[4] < final_time
  }
  return(cond)
}

#' @export
sls_sim.conditioning <- function(
  LL,
  LS,
  cond
)
{

  if (nrow(LS) > 2) {stop('Currently this only works for 1 shift!')}
  if (cond == 2) {stop('Cond 2 is not supported!')}

  for (clade in LS$clade_id)
  {
    L <- LL[[clade]]
    shifts <- sls_sim.get_shifts_info(LS = LS, clade = clade)
    shifts_times <- shifts$when
    shifted_id <- L[(L[, 5] != 0), 3]

    tp <- 0; ts <- shifts_times[1]; tc <- max(LL[[1]][,1])
    if (clade == 1)
    {
      coords_left  <- sign(L[, 3]) == -sign(shifted_id) #subclade does NOT start from here (M1)
      coords_right <- sign(L[, 3]) ==  sign(shifted_id) #subclade does start from here!!! (M2)
      L_left  <- L[coords_left, ]; dim(L_left)  <- c(sum(coords_left) , 5)
      L_right <- L[coords_right,]; dim(L_right) <- c(sum(coords_right), 5)
      colnames(L_left) <- colnames(L_right) <- colnames(L)

      coords_left_cs  <- (L_left[,1]  > shifts_times[1])
      coords_right_cs <- (L_right[,1] > shifts_times[1])
      L_left_cs  <- L_left [coords_left_cs ,]; dim(L_left_cs)  <- c(sum(coords_left_cs), 5) #lineages in M1 born before tshift
      L_right_cs <- L_right[coords_right_cs,]; dim(L_right_cs) <- c(sum(coords_right_cs), 5) #lineages in M2 born before tshift

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
      surv_s        <- sls:::sls_sim.check_survival(
        L = L,
        final_time = tp
      ) #S survives from s to p
    }
  }

  cond1 <- (surv_left_cp && surv_right_cs && surv_s) || (surv_left_cp && surv_right_cp)
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

# UTILITIES ----

#' @export
sls_sim.read_LL <- function(
  LL
)
{
  if (!is.list(LL)) {stop('LL is not a list!!!')}
  LL2 <- LL
  for (i in seq_along(LL))
  {
    L <- LL[[i]]
    LL2[[i]] <- L[L[,3] != 0, ]
  }
  LL2
}

#' @export
sls_sim.get_standard_LS <- function(
  crown_age = 10,
  N0 = 2,
  shift_time = 4
)
{
  testit::assert(crown_age > shift_time)
  LS <- as.data.frame(matrix(0, nrow = 2, ncol = 4))
  LS[, 1] <- c(0, shift_time)
  LS[, 2] <- c(0, 1)
  LS[, 3] <- c(1, 2)
  LS[, 4] <- c(N0, 1)
  LS[1, 1] <- crown_age
  colnames(LS) <- c("birth_time", "motherclade", "clade_id", "N0")

  t0s <- LS[, 1]
  motherclades <- LS[, 2]
  N0s <- LS[, 4]
  testit::assert(
    all(motherclades < seq(from = 1, to = length(motherclades)))
  )
  testit::assert(all(N0s >= 1))
  if (any(N0s[-1] != 1)) {stop('Every subclade should start with 1 species!')}
  testit::assert(all(t0s > 0))

  LS
}

#' @export
sls_sim.get_shifts_info <- function(
  LS = sls:::sls_sim.get_standard_LS(),
  clade
)
{
  when  <- LS[LS[, 2] == clade, 1]
  where <- LS[LS[, 2] == clade, 3]

  info <- data.frame(when = when, where = where)

}

#' @export
sls_sim.get_t0 <- function(
  LS = sls:::sls_sim.get_standard_LS(),
  clade
)
{
  t0 <- LS[LS[, 3] == clade, 1]
}

#' @export
sls_sim.get_N0 <- function(
  LS = sls:::sls_sim.get_standard_LS(),
  clade
)
{
  t0 <- LS[LS[, 3] == clade, 4]
}

#' @export
sls_sim.get_motherclade <- function(
  LS = sls:::sls_sim.get_standard_LS(),
  clade
)
{
  motherclade <- LS[LS[, 3] == clade, 2]
}

#' @export
sls_sim.get_pool <- function(
  L
)
{
  L2 <- sls_sim.cut_L(L)
  right_rows <- L2[, 3] != 0
  if (length(right_rows) == 0) {return(NULL)}
  # L2 <- L[right_rows, ]; dim(L2) <- c(sum(right_rows), ncol(L))
  L3 <- matrix(
    L2[right_rows, ],
    nrow = sum(right_rows),
    ncol = ncol(L2)
  )
  pool <- L3[L3[, 4] == -1, 3]
}

#' @export
sls_sim.cut_L <- function(
  L
)
{
  Nmax <- max(abs(L[, 3]))
  L2 <- L[1:Nmax,]
  dim(L2) <- c(Nmax, ncol(L))
  return(L2)
}
