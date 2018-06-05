dataset = dataset_2
times_matrix <- dataset$times_matrix
lambdas      <- dataset$lambdas
mus          <- dataset$mus
coords       <- times_matrix2t_coordinates(times_matrix = times_matrix)
ti           <- coords$ti
tb           <- coords$tb
ts           <- coords$ts
tf           <- coords$tf

tbranching = tb[1,1]; tshift1 = ts[1,1]; tshift2 = ts[1,2]; N0 <- 1
A <- tshift1 - ti; B <- tbranching - tshift1; C <- tshift2 - tbranching; D <- tf - tshift2
pars <- list(); pars[[1]] <- c(lambdas[1], mus[1]); pars[[2]] <- c(lambdas[2], mus[2]); pars[[3]] <- c(lambdas[3], mus[3])
arrived_at_shift1 <- arrived_at_shift2 <- arrived_at_branching2 <- arrived_at_branching <- shifted1 <- shifted2 <- 0
A <- A*2; B <- B*2; C <- C*2

lineages <- lineages0
lineages <- c(id = 1, N = N0, pars = pars[[1]]);lineages
lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = A);lineages
lineages <- sim_event_shift(lineages = lineages, new_pars = pars[[2]], shifting_id = 1);lineages
lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = B);lineages
lineages <- sim_event_branching(lineages = lineages, branching_id = 1);lineages
lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = C);lineages

lineages_test <- sim_make_branching3(lineages = lineages, branching_lineage = lineages[2,], sim_get_ids(lineages));lineages_test

sim_make_branching3       <- function(lineages, branching_lineage, ids){

  branching_id <- unname(branching_lineage[1])
  branching_N  <- unname(branching_lineage[2])
  lambda       <- unname(branching_lineage[3])
  mu           <- unname(branching_lineage[4])

  new_id_visible   <- (max(ids) + 1)
  new_id_invisible <- (min(ids) - 1)

  new_lineage_visible    <- c(id = new_id_visible  , N = 1, lambda = lambda, mu = mu)
  new_lineage_invisible  <- c(id = new_id_invisible, N = branching_N - 1, lambda = lambda, mu = mu)
  lineages     <- rbind(new_lineage_invisible, lineages, new_lineage_visible)
  lineages[which(lineages[,1]==branching_id), 2] <- 1

  lineages <- sim_tidy_up_lineages(lineages)
  return(lineages)
}
