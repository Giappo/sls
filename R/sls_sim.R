#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
check_survival <- function(L, final_time = 0) {
  if (is.matrix(L))
  {
    cond <- any(L[,4] < final_time)
  }else
  {
    cond <- L[4] < final_time
  }
  return(cond)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
sls_sim <- function(pars, t0, starting_species = c(2,1), cond) {

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

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
sls_sim2 <- function(pars1, age, soc, cond) {

  missnumspec <- c(0,0)
  lambdas <- c(pars1[1], pars1[4])
  mus     <- c(pars1[2], pars1[5])
  Ks      <- c(pars1[3], pars1[6])
  td      <- pars1[7]

  N0s    <- c(soc, 1)
  t0     <- c(age, td)

  testit::assert(all(lambdas >= 0))
  testit::assert(all(mus >= 0))
  testit::assert(soc > 0)
  testit::assert(length(lambdas) == length(mus))
  testit::assert(all(sign(t0) > 0) | all(sign(t0) < 0))
  t0 <- abs(t0)
  Nclades <- length(lambdas)

  shifts <- L.list <- brts.list <- tes.list <- tas.list <- vector("list", Nclades)

  keep_the_sim <- 0
  while (keep_the_sim == 0)
  {#sim loop
    clade <- 0
    clades_born <- c(1, rep(0, Nclades - 1))
    while ((clade <- clade + 1) <= Nclades)
    {#clade loop
      #clade settings
      lambda <- lambdas[clade]
      mu     <- mus[clade]
      N0     <- N0s[clade]
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

    brts.list <- vector("list", Nclades)
    for (clade in 1:Nclades)
    {
      L <- L.list[[clade]]
      if (clade == 1)
      {
        L[L[,4] == td, 4] <- -1
        if (sls:::check_survival(L) & all(L[,3] != 0)) #if L survives
        {
          if (sum(L[,4] == -1) > 1)
          {
            time_points <- unlist(unname(sort(DDD:::L2brts(L, dropextinct = TRUE), decreasing = TRUE)) )
            brts <- -sort(abs(as.numeric(time_points)), decreasing = TRUE)
            if (cond != 0) {brts[1] <- t0[clade]}
            brts.list[[clade]] <- abs(brts)
          }else
          {
            brts <- NULL
          }
        }
      }else
      {
        if (sls:::check_survival(L) & all(L[,3] != 0)) #if L survives
        {
          if (sum(L[,4] == -1) > 1)
          {
            time_points <- unlist(unname(sort(DDD:::L2brts(L, dropextinct = TRUE), decreasing = TRUE)) )
            brts <- -sort(abs(as.numeric(time_points)), decreasing = TRUE)
            brts.list[[clade]] <- abs(brts)
          }
        }
      }
    }; brts.list

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
    LM1cs <- LM1[condLM1cs,]; dim(LM1cs) <- c(sum(condLM1cs), 4) #lineages in M1 born before tshift
    LM2cs <- LM2[condLM2cs,]; dim(LM2cs) <- c(sum(condLM2cs), 4) #lineages in M2 born before tshift

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
    keep_the_sim
  }#sim loop

  out <- list(brts = brts.list, L = L.list)
  return(out)
  # return(L.list)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
DDD_sim <- function (pars, t0, starting_species = c(2,1), cond)
{
  pars1 <- pars
  ddmodel = 1
  age <- t0[1]
  pars <- c(pars1[1], pars1[2], Inf, pars1[3], pars1[4], Inf, t0[2])

  keep_the_sim <- 0
  while (keep_the_sim == 0)
  {
    done = 0
    if (pars[7] > age) {
      stop("The key innovation time is before the crown age of the main clade.")
    }
    if ((pars[1] < pars[2]) | (pars[4] < pars[5])) {
      stop("lambda0 is smaller than mu for one or both clades")
    }
    if (min(pars) < 0) {
      stop("One of the parameters is negative")
    }
    if (!(ddmodel %in% c(1, 1.3, 2, 2.1, 2.2, 2.3, 3, 4, 4.1,
                         4.2))) {
      stop("This diversity-dependence model does not exist or is not implemented")
    }
    while (done == 0) {
      t = rep(0, 1)
      L = matrix(0, 2, 5)
      i = 1
      t[1] = 0
      NM = 2
      NS = 0
      L[1, 1:5] = c(0, 0, -1, -1, 0)
      L[2, 1:5] = c(0, -1, 2, -1, 0)
      linlistM = c(-1, 2)
      linlistS = NULL
      newL = 2
      tinn = age - pars[7]
      ff = DDD:::dd_KI_lamuN(ddmodel, pars, c(NM[i], NS[i]))
      laMN = ff[1]
      muMN = ff[2]
      laSN = ff[3]
      muSN = ff[4]
      denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
      t[i + 1] = t[i] - log(runif(1))/denom
      if (t[i + 1] > tinn & t[i] < tinn) {
        NM[i] = NM[i] - 1
        NS[i] = NS[i] + 1
        linlistS = DDD:::sample2(linlistM, 1)
        L[abs(linlistS), 5] = 1
        linlistM = linlistM[-which(linlistM == linlistS)]
        ff = DDD:::dd_KI_lamuN(ddmodel, pars, c(NM[i], NS[i]))
        laMN = ff[1]
        muMN = ff[2]
        laSN = ff[3]
        muSN = ff[4]
        denom = (laMN + muMN) * NM[i] + (laSN + muSN) *
          NS[i]
        t[i + 1] = tinn - log(runif(1))/denom
      }
      while (t[i + 1] <= age) {
        event = DDD:::sample2(x = 1:4, size = 1, prob = c(laMN * NM[i], muMN * NM[i], laSN * NS[i], muSN * NS[i]))
        i = i + 1
        if (event == 1) {
          ranL = DDD:::sample2(linlistM, 1)
          NM[i] = NM[i - 1] + 1
          NS[i] = NS[i - 1]
          newL = newL + 1
          L = rbind(L, c(t[i], ranL, sign(ranL) * newL, -1, 0))
          linlistM = c(linlistM, sign(ranL) * newL)
        }
        else if (event == 3) {
          ranL = DDD:::sample2(linlistS, 1)
          NM[i] = NM[i - 1]
          NS[i] = NS[i - 1] + 1
          newL = newL + 1
          L = rbind(L, c(t[i], ranL, sign(ranL) * newL, -1, 1))
          linlistS = c(linlistS, sign(ranL) * newL)
        }
        else if (event == 2) {
          ranL = DDD:::sample2(linlistM, 1)
          NM[i] = NM[i - 1] - 1
          NS[i] = NS[i - 1]
          L[abs(ranL), 4] = t[i]
          w = which(linlistM == ranL)
          linlistM = linlistM[-w]
          linlistM = sort(linlistM)
        }
        else if (event == 4) {
          ranL = DDD:::sample2(linlistS, 1)
          NM[i] = NM[i - 1]
          NS[i] = NS[i - 1] - 1
          L[abs(ranL), 4] = t[i]
          w = which(linlistS == ranL)
          linlistS = linlistS[-w]
          linlistS = sort(linlistS)
        }
        if (sum(c(linlistM, linlistS) < 0) == 0 | sum(c(linlistM, linlistS) > 0) == 0) {
          t[i + 1] = Inf
        }
        else {
          ff = DDD:::dd_KI_lamuN(ddmodel, pars, c(NM[i], NS[i]))
          laMN = ff[1]
          muMN = ff[2]
          laSN = ff[3]
          muSN = ff[4]
          denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
          t[i + 1] = t[i] - log(runif(1))/denom
          if (t[i + 1] > tinn & t[i] < tinn) {
            NM[i] = NM[i] - 1
            NS[i] = NS[i] + 1
            ff = DDD:::dd_KI_lamuN(ddmodel, pars, c(NM[i], NS[i]))
            laMN = ff[1]
            muMN = ff[2]
            laSN = ff[3]
            muSN = ff[4]
            linlistS = DDD:::sample2(linlistM, 1)
            L[abs(linlistS), 5] = 1
            linlistM = linlistM[-which(linlistM == linlistS)]
            denom = (laMN + muMN) * NM[i] + (laSN + muSN) *
              NS[i]
            t[i + 1] = tinn - log(runif(1))/denom
          }
        }
      }
      if (sum(c(linlistM, linlistS) < 0) == 0 | sum(c(linlistM,
                                                      linlistS) > 0) == 0) {
        done = 0
      }
      else {
        done = 1
      }
    }
    L[, 1] = age - c(L[, 1])
    notmin1 = which(L[, 4] != -1)
    L[notmin1, 4] = age - c(L[notmin1, 4])
    L[which(L[, 4] == age + 1), 4] = -1
    tes = DDD:::L2phylo(L[, 1:4], dropextinct = T)
    tas = DDD:::L2phylo(L[, 1:4], dropextinct = F)
    tesS = NULL
    tes2 = NULL
    par(mfrow = c(2, 1))
    plot(tes)
    plot(tas)
    cols = c("blue", "red")
    names(cols) = c(0, 1)
    if (length(linlistS) > 0) {
      namesS = paste("t", abs(linlistS), sep = "")
      if (length(linlistS) == 1) {
        m = which(tes$tip.label == namesS)
        b2 = 0
      }
      else if (length(linlistS) > 1) {
        m = ape:::getMRCA(phy = tes, tip = namesS)
        tesS = ape:::extract.clade(phy = tes, node = m)
        b2 = age - ape:::node.depth.edgelength(tes)[m]
      }
      m0 = tes$edge[which(tes$edge[, 2] == m), 1]
      b1 = age - ape:::node.depth.edgelength(tes)[m0]
      tes2 = phytools:::paintSubTree(tes, node = m, state = "1", anc.state = "0",
                                     stem = (pars[7] - b2)/(b1 - b2))
      phytools:::plotSimmap(tes2, cols, lwd = 3, pts = F)
    }
    tasS = NULL
    tas2 = NULL
    allS = which(L[, 5] == 1)
    if (length(allS) > 0) {
      namesS = paste("t", abs(allS), sep = "")
      if (length(allS) == 1) {
        m = which(tas$tip.label == namesS)
        b2 = 0
      }
      else if (length(allS) > 1) {
        m = ape:::getMRCA(phy = tas, tip = namesS)
        tasS = ape:::extract.clade(phy = tas, node = m)
        b2 = age - ape:::node.depth.edgelength(tas)[m]
      }
      m0 = tas$edge[which(tas$edge[, 2] == m), 1]
      b1 = age - ape:::node.depth.edgelength(tas)[m0]
      tas2 = phytools:::paintSubTree(tas, node = m, state = "1", anc.state = "0",
                                     stem = (pars[7] - b2)/(b1 - b2))
      phytools:::plotSimmap(tas2, cols, lwd = 3, pts = F)
    }


    LM <- L[L[,5] == 0, 1:4]; dim(LM) <- c(sum(L[,5] == 0), 4)
    LS <- L[L[,5] == 1, 1:4]; dim(LS) <- c(sum(L[,5] == 1), 4)
    shift_time <- pars[7]

    signS <- unique(sign(LS[,3]))
    condLM1 <- sign(LM[,3]) != signS
    condLM2 <- sign(LM[,3]) == signS
    LM1 <- LM[condLM1,]; dim(LM1) <- c(sum(condLM1), 4)
    LM2 <- LM[condLM2,]; dim(LM2) <- c(sum(condLM2), 4)
    condLM2cs <- (LM2[,1] > shift_time)
    LM2cs <- LM2[condLM2cs,]; dim(LM2cs) <- c(sum(condLM2cs), 4)

    condM1cp <- sls:::check_survival(L = LM1, final_time = 0)
    condM2cp <- sls:::check_survival(L = LM2, final_time = 0)
    condM2cs <- sls:::check_survival(L = LM2cs, final_time = shift_time)
    condS    <- sls:::check_survival(L = LS, final_time = 0)
    testit::assert(condM2cs >= condM2cp)

    cond1 <- (condM1cp && condM2cs && condS) || (condM1cp && condM2cp)
    cond2 <- (condM1cp && condM2cp && condS)
    cond3 <- (condM1cp && condM2cs && condS)

    #conditioning
    keep_the_sim <- (cond == 0) + (cond == 1) * cond1 + (cond == 2) * cond2 + (cond == 3) * cond3; keep_the_sim
  }
  out = list(tes = tes, tas = tas, L = L, tesS = tesS, tasS = tasS,
             tes2 = tes2, tas2 = tas2)
  return(out)
}


