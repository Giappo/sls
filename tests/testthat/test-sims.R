context("simulations")

test_that( "simulated trees respect conditional instances", {

  testthat::skip('This test includes sls_sim, which I want to replace.')

  cond <- 1
  age <- 8
  N0 <- 2

  simpars <- c(0.5, 0.4, 0.2, 0.1)
  t_d <- 4
  pars1 <- c(simpars[1], simpars[2], Inf, simpars[3], simpars[4], Inf, t_d)

  maxseed <- 3 + (7 * ribir::is_on_travis()); s <- 1; set.seed(s)
  for (s in 1:maxseed)
  {
    print(s)
    set.seed(s)
    for (cond in 1:3)
    {
      sim <- sls::sls_sim(pars1 = pars1, age = age, N0 = N0, cond = cond)
      L <- sim$L
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
      sim <- sls::sls_sim(pars1 = pars1, age = age, N0 = N0, cond = cond)
      brts <- sim$brts
      L1 <- sim$L[[1]]
      alive_at_the_present <- (L1[,4] < 0)
      left  <- sum(sign(L1[alive_at_the_present, 3]) ==  1)
      right <- sum(sign(L1[alive_at_the_present, 3]) == -1)
      both_crown_present <- sign(left * right)
      testthat::expect_true( (both_crown_present * brts[[1]][1]) == (both_crown_present * age) )
    }
  }
})
