# RJCB: Well written!

#basic sim function
#' @export
sim_bd                   <- function(pars, time, N0 = 1){ #<- function(lambda, mu, ti, tf, N0 = 1){

  ti <- 0; tf <- time;
  lambda = pars[1]; mu = pars[2]
  N <- N0
  t <- ti #continous time
  while(t < tf && N > 0)
  {
    total_rate <- N * (lambda + mu)
    deltaT <- rexp(n = 1, rate = total_rate)
    if ((t + deltaT) < tf)
    {
      deltaN <- DDD:::sample2(c(-1,1), size = 1, replace = F, prob = c(mu,lambda))
      N <- N + deltaN
      t <- t + deltaT
    }else
    {
      t <- tf
    }
  }
  return(N)
}

#sim modules lvl 1
#' @export
sim_produce_bad_lineages <- function(lineages){
  #define bad_lineages aka the output to return if something goes wrong
  bad_lineages <- lineages
  if (is.matrix(lineages)==0) #if there is only one
  {
    bad_lineages[2] <- 0
  }else #if there is more
  {
    bad_lineages[,2] <- 0
  }
  return(bad_lineages)
}
#' @export
sim_get_ids              <- function(lineages){
  #get lineages ids
  if (!is.matrix(lineages)) #if there is only one
  {
    ids <- lineages[1]
  }else #if there is more
  {
    ids <- lineages[,1]
  }
  return(ids)
}
#' @export
sim_get_changing_lineage <- function(lineages, ids, changing_id){
  #determine which lineage is going to branch
  if (!is.matrix(lineages)) #if there is only one
  {
    branching_lineage <- lineages
    if (changing_id != ids)
    {
      branching_lineage[2] = 0
    }
  }else #if there is more
  {
    if (!(changing_id %in% ids))
    {
      branching_lineage <- lineages[1,]
      branching_lineage[2] = 0
    }else
    {
      corresponding_id  <- which(ids==changing_id)
      branching_lineage <- lineages[corresponding_id,]
    }
  }
  return(branching_lineage)
}
#' @export
sim_make_branching       <- function(lineages, branching_lineage, ids){

  lambda <- branching_lineage[3]; mu <- branching_lineage[4]
  new_id <- (max(ids) + 1)
  ids <- c(ids,new_id)
  new_lineage <- c(id = new_id, N = 1, lambda = lambda, mu = mu)

  lineages <- rbind(lineages, new_lineage)
  lineages <- sim_tidy_up_lineages(lineages)
  return(lineages)
}
#' @export
sim_make_shifting        <- function(lineages, shifting_lineage, ids, new_pars){

  shifting_id <- shifting_lineage[1]
  old_pars <- c(shifting_lineage[3], shifting_lineage[4])
  new_pars <- new_pars
  old_id <- shifting_id
  new_id <- (min(ids) - 1)
  ids <- c(ids, new_id)
  shifted_lineage <- c(id = new_id, N = (shifting_lineage[2] - 1), lambda = old_pars[1], mu = old_pars[2])
  shifting_lineage<- c(id = old_id, N = 1                        , lambda = new_pars[1], mu = new_pars[2])
  if (is.matrix(lineages)==0) #if there is only one
  {
    lineages <- rbind(shifting_lineage, shifted_lineage)
    colnames(lineages) <- c("id","N","lambda","mu")
  }else #if there is more
  {
    lineages[which(lineages[,1]==old_id),] <- shifting_lineage
    lineages <- rbind(lineages, shifted_lineage)
    colnames(lineages) <- c("id","N","lambda","mu")
  }

  lineages <- sim_tidy_up_lineages(lineages)
  return(lineages)
}

#sim modules lvl 2
#' @export
sim_event_branching      <- function(lineages, branching_id){

  #define bad_lineages aka the output to return if something goes wrong
  bad_lineages <- sim_produce_bad_lineages(lineages = lineages)

  #get lineages ids
  ids <- sim_get_ids(lineages = lineages)

  #determine which lineage is going to branch (if branching_id is not present in ids, branching lineage will have N = 0)
  branching_lineage <- sim_get_changing_lineage(lineages = lineages, ids = ids, changing_id = branching_id)

  #if that lineage is not present, return lineages without any modification (will lead to a 0 ok score)
  if (is.na(branching_lineage[2])) {return(bad_lineages)}

  #if there are no individuals, return lineages without any modification (will lead to a 0 ok score)
  if (branching_lineage[2] == 0) {return(bad_lineages)}

  #if the designated lineages is present and it has at least one individual then BRANCH!
  lineages <- sim_make_branching(lineages = lineages, ids = ids, branching_lineage = branching_lineage)

  return(lineages)
}
#' @export
sim_event_shift          <- function(lineages, shifting_id, new_pars){

  #define bad_lineages aka the output to return if something goes wrong
  bad_lineages <- sim_produce_bad_lineages(lineages = lineages)

  #get lineages ids
  ids <- sim_get_ids(lineages = lineages)

  #determine which lineage is going to shift (if shifting_id is not present in ids, shifting lineage will have N = 0)
  shifting_lineage <- sim_get_changing_lineage(lineages = lineages, ids = ids, changing_id = shifting_id)

  #if that lineage is not present, return lineages without any modification (will lead to a 0 ok score)
  if (is.na(shifting_lineage[2])) {return(bad_lineages)}

  #check if the shifting lineage has N=0 (it should never happen, but better safe than sorry)
  if (shifting_lineage[2] == 0) {return(bad_lineages)}

  #if it is all ok, update the lineages matrix to account for the shift effects
  lineages <- sim_make_shifting(lineages = lineages, shifting_lineage = shifting_lineage, ids = ids, new_pars = new_pars)

  return(lineages)
}
#' @export
sim_evolve_lineages      <- function(lineages, time_before_next_event){

  if (is.matrix(lineages)==0)
  { #if there is only one
    ids <- lineages[1]
    lineages[2] <- sim_bd(pars = c(lineages[3], lineages[4]) , time = time_before_next_event, N0 = lineages[2] )
  }else
  { #if there is more
    ids      <- lineages[,1]
    Ns       <- lineages[,2]
    lambdas  <- lineages[,3]
    mus      <- lineages[,4]
    for (i in ids)
    {
      ii <- which(ids==i)
      New_N <- sim_bd(pars = c(lambdas[ii], mus[ii]) , time = time_before_next_event, N0 = Ns[ii] )
      lineages[ii,2]  <- New_N
    }
  }

  lineages <- sim_tidy_up_lineages(lineages)
  return(lineages)
}
#' @export
sim_check_ok_condition   <- function(lineages, Ntips){
  ok_dead <- 0; ok_survivors <- 0; ok_correct_number_of_tips <- 0

  if (is.matrix(lineages))
  {
    if (length(lineages[,1] >  0) > 0)
    {
      survivors    <- which(lineages[,1] >  0)
      ok_survivors <- prod( lineages[survivors, 2] == 1 )
      ok_correct_number_of_tips = (length(survivors) == Ntips)
    }
    if (length(lineages[,1] <= 0) > 0)
    {
      deads        <- which(lineages[,1] <= 0)
      ok_dead      <- prod( lineages[deads, 2]     <= 0 )
    }
  }else
  { #if there is only one lineage
    ok_survivors <- (lineages[2] == 1)
    ok_dead <- 1
  }

  ok <- ok_survivors * ok_dead * ok_correct_number_of_tips
  return(ok)
}
#' @export
sim_tidy_up_lineages     <- function(lineages){

  if (is.matrix(lineages))
  {
    tidy_lineages <- lineages[order(lineages[,1],decreasing = T),]

    visibles   <- tidy_lineages[,1]> 0
    invisibles <- tidy_lineages[,1]<=0
    if ((sum(visibles) + sum(invisibles)) > 1)
    {
      rownames(tidy_lineages)[visibles]   <- rep("visible"  , sum(visibles))
      rownames(tidy_lineages)[invisibles] <- rep("invisible", sum(invisibles))
    }
    lineages <- tidy_lineages
  }
  return(lineages)
}

#main sim
#' @export
sim_custom               <- function(lambdas, mus, ti, tb, ts, tf, N0 = 1, input_check = TRUE){

  # Checks if the input data is coherent
  if (input_check == TRUE)
  {
    testit::assert(length(lambdas) == length(mus))
    testit::assert(tf > ti)
    testit::assert(N0 > 0)
    testit::assert(nrow(tb) == 2 || is.null(nrow(tb)))
    testit::assert(nrow(ts) == 2 || is.null(nrow(ts)))
    coherent_input <- check_input_data_coherence(ti = ti, tf = tf, tb = tb, ts = ts)
    if (coherent_input == 0){stop("Input data are incoherent")}
  }

  # Combine the event matrix "times_matrix": first line are time points, second line are ids. Negative ids are shifts.
  times_matrix  <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
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
      lineages <- sim_event_shift(lineages = lineages, new_pars = pars[[regime]], shifting_id = abs(times_matrix[2,t]))#; lineages
    }
    if (branch)
    {
      lineages <- sim_event_branching(lineages = lineages, branching_id = abs(times_matrix[2,t]))
    }
  }

  ok <- sim_check_ok_condition(lineages = lineages, Ntips = Ntips); ok
  return(list(ok = ok, lineages = lineages))
}
