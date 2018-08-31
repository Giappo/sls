#' @inheritParams default_params_doc
#' @export
HH <- function(lambda, mu, t) {
  LL  <- exp((mu - lambda) * t)
  out <- ((lambda - mu)^2 * LL) / (lambda - mu * LL)^2
  return(out)
}

#' @inheritParams default_params_doc
#' @export
H_lik <- function(brts, shifts, lambdas, mus) {

  lambda1 <- 0.4; mu1 <- 0.1
  lambda2 <- 0.3; mu2 <- 0.05
  brts <- c(-10, -5, -3)
  shifts <- c(-4)

  lambda1 <- lambdas[1]; lambda2 <- lambdas[2]
  mu1 <- mus[1]; mu2 <- mus[2]

  age <- max(abs(brts))
  brts <- brts + age
  shifts <- shifts + age#; shifts <- -shifts
  times <- sort(c(brts, shifts))
  times <- times * (-((times %in% shifts)-0.5)*2)
  N0 <- 2
  kvec <- N0 + c(cumsum(sign(times)))

  Ht <- rep(0, length(times) + 1)
  Ht[1] <- Ht[2] <- HH(lambda = lambda1, mu = mu1, t = age - times[1])
  for (t in 2:length(times))
  {
    Ht[t + 1] <- HH(lambda = lambda1, mu = mu1, t = age - abs(times[t])) ^ sign(times[t])
  }; Ht
  L1 <- prod(Ht) * prod(kvec) * lambda^sum(sign(times[-1]) == 1)
  L2 <- prod(HH(lambda = lambda2, mu = mu2, t = age - abs(shifts)))
  lik <- L1 * L2
}

#' @inheritParams default_params_doc
#' @export
sim_get_changing_lineage_random <- function(lineages, ids, changing_id){
  #determine which lineage is going to branch
  if (!is.matrix(lineages)) #if there is only one
  {
    branching_lineage <- lineages
  }else #if there is more
  {
    Ns <- lineages[,2]
    ids2 <- ids[ids > 0 & Ns > 0]
    Ns2  <- Ns[ids > 0 & Ns > 0]
    if (length(Ns2) == 1)
    {
      corresponding_id  <- which(ids == ids2)
      branching_lineage <- lineages[corresponding_id,]
    }else if (length(Ns2) > 1)
    {
      corresponding_id  <- which(ids == sample(x = ids2, size = 1, prob = Ns2))
      branching_lineage <- lineages[corresponding_id,]
    }else
    {
      branching_lineage <- lineages[1,]
    }
  }
  return(branching_lineage)
}

#' @inheritParams default_params_doc
#' @export
sim_custom_custom            <- function(dataset, N0 = 1, input_check = TRUE,
                                         branching_function = sls::sim_event_branching,
                                         shift_function = sls::sim_event_shift,
                                         branching_lineage_selecter = sls::sim_get_changing_lineage,
                                         branching_maker = sls::sim_make_branching,
                                         shifting_lineage_selecter = sls::sim_get_changing_lineage,
                                         shifting_maker = sls::sim_make_shifting){

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus
  coords       <- sls::times_matrix2t_coordinates(times_matrix = times_matrix)
  ti           <- coords$ti
  tb           <- coords$tb
  ts           <- coords$ts
  tf           <- coords$tf

  # Checks if the input data is coherent
  if (input_check == TRUE)
  {
    testit::assert(length(lambdas) == length(mus))
    testit::assert(tf > ti)
    testit::assert(N0 > 0)
    testit::assert(nrow(tb) == 2 || is.null(nrow(tb)))
    testit::assert(nrow(ts) == 2 || is.null(nrow(ts)))
    coherent_input <- check_input_data_coherence(dataset = dataset, N0 = N0)
    if (coherent_input == 0){stop("Input data are incoherent")}
  }

  # Combine the event matrix "times_matrix": first line are time points, second line are ids. Negative ids are shifts.
  pars <- list(); for (i in 1:length(lambdas)) {pars[[i]] <- c(lambdas[i], mus[i])}

  nbranches <- 0; if(!is.null(ncol(tb))){nbranches <- ncol(tb)};
  Ntips <-  N0 + nbranches; r <- rep(regime <- 1, Ntips); n <- 1; t <- 2; shift <- branch <- 0
  lineages <- c(id = 1, N = N0, pars = pars[[1]])
  for (t in 2:length(times_matrix[1,]))
  {
    time_before_next_event <- times_matrix[1,t] - times_matrix[1,t-1]
    lineages <- sim_evolve_lineages(lineages = lineages, time_before_next_event = time_before_next_event);

    #the sign of times_matrix's second line determines if it is a shift (-1), a branching (+1), or neither of the two (0)
    shift    <- (sign(times_matrix[2,t]) < 0)
    branch   <- (sign(times_matrix[2,t]) > 0)
    if (sim_check_id_presence(lineages = lineages, id = abs(times_matrix[2,t])) == 0)
    {
      return(list(ok = 0, lineages = lineages))
    }
    if (shift)
    {
      regime   <- regime + 1
      lineages <- shift_function(lineages = lineages,
                                 new_pars = pars[[regime]],
                                 shifting_id = abs(times_matrix[2,t]),
                                 shifting_lineage_selecter = shifting_lineage_selecter,
                                 shifting_maker = shifting_maker)
    }
    if (branch)
    {
      lineages <- branching_function(lineages = lineages,
                                     branching_id = abs(times_matrix[2,t]),
                                     branching_lineage_selecter = branching_lineage_selecter,
                                     branching_maker = branching_maker)
    }
  }

  ok <- sim_check_ok_condition(lineages = lineages, Ntips = Ntips); ok
  return(list(ok = ok, lineages = lineages))
}
