context("simulations")

test_that( "simulated trees respect conditional instances", {

  lambdas <- c(0.5, 0.4)
  mus     <- c(0.2, 0.1)
  pars    <- c(lambdas[1], mus[1], lambdas[2], mus[2])
  t0      <- c(10, 6)
  starting_species <- c(2,1)
  maxseed <- 10; s <- cond <- 1; set.seed(s)
  for (s in 1:maxseed)
  {
    print(s)
    set.seed(s)
    for (cond in 1:3)
    {
      L <- sls::sls_sim(pars, t0 = t0, starting_species = starting_species, cond = cond)$L
      checkM <- unique(sign(L[[1]][(L[[1]][,4] == -1),3])); checkM
      checkS <- unique(sign(L[[2]][(L[[2]][,4] == -1),3])); checkS

      cond1 <- (length(unique(c(checkM, checkS))) == 2) #among M and S there are descendants of both crown species
      cond2 <- ((length(checkM) + length(checkS)) == 3) #both crown species have descendants at the present in the main clade. Subclade survives
      cond3 <- (length(checkS) > 0) && (length(unique(c(checkM, checkS))) == 2) #among M and S there are descendants of both crown species. Subclade survives

      ok <- (cond == 0) +
            (cond == 1) * cond1 +
            (cond == 2) * cond2 +
            (cond == 3) * cond3
      testthat::expect_true(ok == 1)
    }
  }

  for (s in (maxseed + 1):(2 * maxseed))
  {
    set.seed(s)
    for (cond in 1:3)
    {
      test <- sls::sls_sim(pars, t0 = t0, starting_species = starting_species, cond = cond);test$L[[2]]; test$brts[[2]]
      brts <- test$brts
      L1 <- test$L[[1]]
      alive_at_the_present <- (L1[,4] < 0)
      left  <- sum(sign(L1[alive_at_the_present, 3]) ==  1)
      right <- sum(sign(L1[alive_at_the_present, 3]) == -1)
      both_crown_present <- sign(left * right)
      testthat::expect_true((both_crown_present * brts[[1]][1]) == -(both_crown_present * t0[1]))
      testthat::expect_true(brts[[2]][1] == -t0[2])
    }
  }
})
