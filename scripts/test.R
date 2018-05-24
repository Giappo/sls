# source("R/functions_likelihood.R"); source("R/functions_sim.R"); source("R/functions_utility.R")
# lineages = c(id = 1, N = 1, lambda = 0.2, mu = 0); lineages
# ####
# branching_id = 1
# lineages <- sim_event_branching(lineages = lineages, branching_id = branching_id);lineages
# ####
# shifting_id = 2; new_pars = c(lambda = 0.8, mu = 0.4)
# lineages <- sim_event_shift(lineages = lineages, shifting_id = shifting_id, new_pars = new_pars);lineages
# ####
# time_before_next_event <- 1
# lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = time_before_next_event);lineages
#
# sim_check_ok_condition(lineages, Ntips = 1)
#
# ###
#
# test_this <- function(){
# time_before_next_event <- 1
#
# lineages = c(id = 1, N = 1, lambda = 0.2, mu = 0); lineages
# lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = time_before_next_event);lineages
# branching_id = 1
# lineages <- sim_event_branching(lineages = lineages, branching_id = branching_id);lineages
# lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = time_before_next_event);lineages
# shifting_id = 2; new_pars = c(lambda = 0.1, mu = 0)
# lineages <- sim_event_shift(lineages = lineages, shifting_id = shifting_id, new_pars = new_pars);lineages
# lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = time_before_next_event);lineages
# branching_id = 1
# lineages <- sim_event_branching(lineages = lineages, branching_id = branching_id);lineages
# lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = time_before_next_event);lineages
# branching_id = 1
# lineages <- sim_event_branching(lineages = lineages, branching_id = branching_id);lineages
# lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = time_before_next_event);lineages
# ok <- sim_check_ok_condition(lineages = lineages, Ntips = 4)}
# pippo <- test_this(); pippo
#
# oks <- rep(NA,maxN <- 1000000)
# for (i in 1:maxN)
# {
#   oks[i] <- test_this()
# };print(out <- sum(oks)/maxN)
# lik_custom(lambdas = c(0.2, 0.1), mus = c(0,0),ti = 0,)
# prod(exp(-6*c(0.2, 0.1)))
