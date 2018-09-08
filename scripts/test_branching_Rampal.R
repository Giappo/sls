rm(list = ls())
library(sls)
dataset = dataset_2
times_matrix <- dataset$times_matrix
lambdas      <- dataset$lambdas
mus          <- dataset$mus

pars <- list(); pars[[1]] <- c(lambdas[1], mus[1]); pars[[2]] <- c(lambdas[2], mus[2]); pars[[3]] <- c(lambdas[3], mus[3])
Nsims <- 1E6; A <- 3; B <- 4; C <-5; N0 <- 1
# #test 1
# oks <- rep(NA, Nsims)
# for (i in 1:Nsims){
#   lineages  <- c(id = 1, N = N0, pars = pars[[1]]);lineages
#   lineages  <- sim_evolve_lineages(lineages = lineages, time_before_next_event = A);lineages
#   lineages  <- sim_event_shift(lineages = lineages, shifting_id = 1, new_pars = pars[[1]]);lineages
#   lineages  <- sim_evolve_lineages(lineages = lineages, time_before_next_event = B);lineages
#   oks[i]    <- sim_check_ok_condition(lineages = lineages, Ntips = 1)
# };sim <- sum(oks)/Nsims; sim
#
# lik <- lik_Pi(lambda = pars[[1]][1], mu = pars[[1]][2], time = A) *
#        lik_Pi(lambda = pars[[1]][1], mu = pars[[1]][2], time = B) *
#        lik_single_event(lambda = pars[[1]][1], mu = pars[[1]][2], time1 = A, time2 = B);lik
# std <- get_std2(oks = oks, lik_result = lik, sim_result = sim)$std_max
# test_result1 <- list(lik = lik, sim = sim, std = std);test_result1
# lapply(test_result1[1:3], write, "basic_tests.txt", append = TRUE)

#test 2
oks <- rep(NA, Nsims)
for (i in 1:Nsims){
  lineages <- c(id = 1, N = N0, pars = pars[[1]]);lineages
  lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = A);lineages
  if (lineages[2] == 0)
  {
    oks[i] <- 0
  }
  else
  {
    lineages  <- sim_event_branching(lineages = lineages, branching_id = 1);lineages
    lineages  <- lineages[which(lineages[,1] == 1),]; lineages
    lineages  <- sim_evolve_lineages(lineages = lineages, time_before_next_event = B);lineages
    oks[i]    <- sim_check_ok_condition(lineages = lineages, Ntips = 1)
  }
};sim <- sum(oks)/Nsims; sim

lik <- unname(lik_Pi(lambda = pars[[1]][1], mu = pars[[1]][2], time = A+B));lik
std <- get_std2(oks = oks, lik_result = lik, sim_result = sim)$std_max
test_result2 <- list(lik = lik, sim = sim, std = std);test_result2
lapply(test_result2[1:3], write, "basic_tests.txt", append = TRUE)

#test 3
oks <- rep(NA, Nsims); i <- 1
for (i in 1:Nsims){
  lineages <- c(id = 1, N = N0, pars = pars[[1]]);lineages
  lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = A);lineages
  if (lineages[2] == 0){oks[i] <- 0}
  else
  {
    lineages  <- sim_event_branching(lineages = lineages, branching_id = 1);lineages
    lineages  <- lineages[which(lineages[,1] == 1),]; lineages
    lineages  <- sim_evolve_lineages(lineages = lineages, time_before_next_event = B);lineages
    if (lineages[2] == 0){oks[i] <- 0}
    else
    {
      lineages  <- sim_event_shift(lineages = lineages, shifting_id = 1, new_pars = pars[[1]]);lineages
      lineages  <- sim_evolve_lineages(lineages = lineages, time_before_next_event = C);lineages
      oks[i]    <- sim_check_ok_condition(lineages = lineages, Ntips = 1)
    }
  }
};sim <- sum(oks)/Nsims; sim

lik <- lik_Pi(lambda = pars[[1]][1], mu = pars[[1]][2], time = A+B) *
  lik_Pi(lambda = pars[[1]][1], mu = pars[[1]][2], time = C) *
  lik_single_event(lambda = pars[[1]][1], mu = pars[[1]][2], time1 = B, time2 = C);lik
lik2 <- lik_Pi(lambda = pars[[1]][1], mu = pars[[1]][2], time = A) *
  lik_Pi(lambda = pars[[1]][1], mu = pars[[1]][2], time = B) *
  lik_single_event(lambda = pars[[1]][1], mu = pars[[1]][2], time1 = A, time2 = B + C) *
  lik_Pi(lambda = pars[[1]][1], mu = pars[[1]][2], time = C) *
  lik_single_event(lambda = pars[[1]][1], mu = pars[[1]][2], time1 = B, time2 = C);lik2
std <- get_std2(oks = oks, lik_result = lik, sim_result = sim)$std_max
test_result3 <- list(lik = lik, sim = sim, std = std);test_result3
lapply(test_result3[1:3], write, "basic_tests.txt", append = TRUE)



