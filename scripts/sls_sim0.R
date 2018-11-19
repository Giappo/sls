#' Simulation code (old)
#' @inheritParams default_params_doc
#' @return result
#' @export
sls_sim0 <- function(pars, t0, starting_species = c(2,1), cond) {

  lambdas <- c(pars[1], pars[3])
  mus     <- c(pars[2], pars[4])

  testit::assert(all(lambdas >= 0))
  testit::assert(all(mus >= 0))
  testit::assert(all(starting_species > 0))
  testit::assert(length(lambdas) == length(mus))
  testit::assert(length(lambdas) == length(starting_species))
  testit::assert(all(sign(t0) > 0) | all(sign(t0) < 0))
  t0 <- abs(t0)
  Nclades <- length(lambdas)

  shifts <- L.list <- brts.list <- tes.list <- tas.list <- vector("list", Nclades)

  keep_the_sim <- 0
  while (keep_the_sim == 0)
  {
    clade <- 0
    clades_born <- c(1, rep(0, Nclades - 1))
    while ((clade <- clade + 1) <= Nclades)
    {#clade loop
      #clade settings
      lambda <- lambdas[clade]
      mu     <- mus[clade]
      N0     <- starting_species[clade]
      age    <- t0[clade]
      tshift <- t0[t0 < age];
      tshift <- c(tshift, -2)
      shifts[[clade]] <- tshift
      shifts_already_occurred <- 0

      #initial time input
      total_count <- N0
      pool <- 1:N0
      N <- N0
      t <- age
      L <- matrix(0, nrow = 1e6, 4)
      L[,4] <- -1
      L[,3] <- 0
      L[1,1:4] <- c(t, 0,  1, -1)
      if (N0 == 2){L[2,1:4] <- c(t, 1, -2, -1)}; head(L)
      if (clade == 1) {shifted <- 0} else if (clade == 2) {L[1,3] <- sign(shifted)}
      pool <- L[,3][1:N0]

      while (t > 0)
      {#time loop
        # print(pool)
        N <- length(pool)
        total_rate <- N * (lambda + mu)
        if (total_rate > 0)
        {
          deltaT <- rexp(1, rate = total_rate)
          deltaN <- sample(c(-1, 1), size = 1, prob = c(mu, lambda))
          t <- t - deltaT
          if (t < tshift[shifts_already_occurred + 1] & t > 0)
          {
            t <- tshift[shifts_already_occurred + 1]
            if (N > 1) {shifted <- sample(pool, replace = FALSE, size = 1)}else{shifted <- pool}
            L[abs(shifted), 4] <- t
            pool <- pool[pool != shifted]
            shifts_already_occurred <- shifts_already_occurred + 1
          }else
          {
            if (deltaN > 0 & t > 0)
            {
              if (N > 1) {parents <- sample(pool, replace = FALSE, size = deltaN)}else{parents <- pool}
              new_interval <- (total_count + 1):(total_count + deltaN)
              L[new_interval, 1] <- t
              L[new_interval, 2] <- parents
              L[new_interval, 3] <- abs(new_interval) * sign(parents)
              pool <- c(pool, abs(new_interval) * sign(parents) )
              total_count <- total_count + deltaN
            }
            if (deltaN < 0 & t > 0)
            {
              if (N > 1) {dead <- sample(pool, replace = FALSE, size = 1)}else{dead <- pool}
              L[abs(dead), 4] <- t
              pool <- pool[pool != dead]
            }
          }
        }else
        {
          t <- 0
        }
      }#time loop

      #present time output
      t <- 0
      L <- L[(1:total_count),]
      dim(L) <- c(total_count, 4)

      #clade output
      L.list[[clade]] <- L
    }; L.list#clade loop

    if (cond) {L.list[[1]][1,1] <- t0[1]}

    # par(c(2,1))
    for (clade in 1:Nclades)
    {
      # brts.list[[clade]] <- NULL
      L <- L.list[[clade]]
      if (sls:::check_survival(L) & all(L[,3] != 0))
      {
        if (sum(L[,4] < 0) == 1)
        {
          brts <- -t0[clade]
        }else
        {
          time_points <- unlist(unname(sort(DDD:::L2brts(L, dropextinct = TRUE), decreasing = TRUE)) )
          brts0 <- -sort(abs(as.numeric(time_points)), decreasing = TRUE)
          if (starting_species[clade] == 1) {brts <- c(-age, brts0)} else {brts <- brts0}
          # plot(tes <- DDD:::L2phylo(L, dropextinct = TRUE))
        }
        # plot(tas <- DDD:::L2phylo(L, dropextinct = FALSE))
        brts.list[[clade]] <- brts
      }
    }

    LM <- L.list[[1]]
    shift_time <- shifts[[1]][1]
    shifted <- (LM[,4] == shift_time)
    shifted_id <- LM[shifted, 3]
    condLM1 <- sign(LM[,3]) != sign(shifted_id) #subclade doesn't start from here (M1)
    condLM2 <- sign(LM[,3]) == sign(shifted_id) #subclade does start from here!!! (M2)
    LM1 <- LM[condLM1,]; dim(LM1) <- c(sum(condLM1), 4)
    LM2 <- LM[condLM2,]; dim(LM2) <- c(sum(condLM2), 4)
    condLM1cs <- (LM1[,1] > shift_time)
    condLM2cs <- (LM2[,1] > shift_time)
    LM1cs <- LM1[condLM1cs,]; dim(LM1cs) <- c(sum(condLM1cs), 4)
    LM2cs <- LM2[condLM2cs,]; dim(LM2cs) <- c(sum(condLM2cs), 4)

    LS <- L.list[[2]]
    condM1cp <- sls:::check_survival(L = LM1, final_time = 0)
    condM2cp <- sls:::check_survival(L = LM2, final_time = 0)
    condM1cs <- sls:::check_survival(L = LM1cs, final_time = shift_time)
    condM2cs <- sls:::check_survival(L = LM2cs, final_time = shift_time)
    condS    <- sls:::check_survival(L = LS, final_time = 0)
    testit::assert(condM1cs >= condM1cp)
    testit::assert(condM2cs >= condM2cp)

    cond1 <- (condM1cp && condM2cs && condS) || (condM1cp && condM2cp)
    cond2 <- (condM1cp && condM2cp && condS)
    cond3 <- (condM1cp && condM2cs && condS)

    #conditioning
    keep_the_sim <- (cond == 0) + (cond == 1) * cond1 + (cond == 2) * cond2 + (cond == 3) * cond3
    # keep_the_sim <- keep_the_sim * condM1cs #we require that M1 arrives at the shift
  }

  # out <- list(brts = brts.list, tes = tes.list, tas = tas.list, L = L.list)
  # out <- list(L = L.list)
  out <- list(brts = brts.list, L = L.list)
  return(out)
  # return(L.list)
}
