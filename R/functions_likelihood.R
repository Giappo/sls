#BASIC FUNCTIONS

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_pt  <- function (lambda, mu, time){
  Lambda <- exp((mu-lambda)*(time))
  out    <- (lambda==mu)*( 1/(1+lambda*(time)) ) +
    (lambda!=mu)*( (lambda-mu)/(lambda-mu*Lambda*(lambda!=mu)) )
  return(out)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_ut  <- function (lambda, mu, time){
  Lambda <- exp((mu-lambda)*(time))
  out    <- (lambda==mu)*( lambda*(time)/(1+lambda*(time)) ) +
    (lambda!=mu)*( (lambda - lambda*Lambda) / (lambda - mu*Lambda*(lambda!=mu)) )
  return(out)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_Pi  <- function (lambda, mu, time){
  out <- lik_pt(lambda = lambda, mu = mu, time = time)*(1-lik_ut(lambda = lambda, mu = mu, time = time))
  return(out)
}

#LIKELIHOODS

#' Basic denominator likelihood module. If squared yields branching contribution, otherwise is shift contribution.
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_single_event <- function(lambda, mu, time1, time2){ #if squared describes the linkage between two consecutives full lines, if its exponent is one it describes instead a single-lineage shift.

  lik <- 1 - lik_ut(lambda = lambda, mu = mu, time = time1) *
    (1 - lik_pt(lambda = lambda, mu = mu, time = time2))

  lik <- lik^(-1)

  return(lik)
}

#' Custom likelihood function.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom      <- function(dataset, N0 = 1, input_check = 1){

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus
  coords       <- times_matrix2t_coordinates(times_matrix = times_matrix)
  ti           <- coords$ti
  tb           <- coords$tb
  ts           <- coords$ts
  tf           <- coords$tf
  Ntimepoints  <- ncol(times_matrix)

  if (input_check == 1){
    coherent_input <- check_input_data_coherence(dataset = dataset, N0 = N0)
    if (coherent_input == 0){stop("Input data are incoherent")}
  }

  if (!is.null(ncol(tb))){ rownames(tb) <- c("time", "who") }
  if (!is.null(ncol(ts))){ rownames(ts) <- c("time", "who") }
  nbranches <- 0; if(!is.null(ncol(tb))){nbranches <- ncol(tb)};
  Ntips <-  N0 + nbranches
  fathers <- unname( c(rep(0, N0), tb[2,]) ); sons <- 1:Ntips
  lik_den <- 1; r <- rep(regime <- 1, Ntips); n <- 1; t <- 2; shift <- branch <- 0
  for (t in 2:(Ntimepoints - 1))
  {
    shift  <- (sign(times_matrix[2,t]) < 0)
    branch <- (sign(times_matrix[2,t]) > 0)

    lineages_affected_by_the_event <- rep(0,n); lineages_affected_by_the_event[abs(times_matrix[2,t])] <- 1
    expo <- unname(2 - lineages_affected_by_the_event * shift)
    # print(lineages_affected_by_the_event)

    time_interval1 <- times_matrix[1,t] - times_matrix[1,t-1]; time_interval2 <- tf - times_matrix[1,t]
    per_lineage_outcome <- lik_single_event(time1 = time_interval1, time2 = time_interval2, lambda = lambdas[r[1:n]], mu = mus[r[1:n]])^expo
    lik_den <- lik_den * prod(per_lineage_outcome)

    if(shift)
    {
      updated <- update_regimes(t = t, regime = regime, times_matrix = times_matrix, r = r, n = n)
      r <- updated$r
      regime <- updated$regime
    }
    n <- n + branch
    # t <- t + 1; t #remember to comment this
    # print(r); #remember to comment this
  }#; r

  lik_nom <- 1; n <- 1; t <- 2
  r <- rep(1, Ntips)
  shift <- branch <- 0
  regime <- 1
  for (t in 2:Ntimepoints)
  {
    shift  <- (sign(times_matrix[2,t]) < 0)
    branch <- (sign(times_matrix[2,t]) > 0)

    # print(r[1:n])
    # print(lambdas[r[1:n]])
    time_interval <- times_matrix[1,t] - times_matrix[1,t-1]
    per_lineage_outcome <- lik_Pi(lambda = lambdas[r[1:n]], mu = mus[r[1:n]], time = time_interval)
    lik_nom <- lik_nom * prod(per_lineage_outcome)

    if(shift)
    {
      updated <- update_regimes(t = t, regime = regime, times_matrix = times_matrix, r = r, n = n)
      r <- updated$r
      regime <- updated$regime
    }
    n <- n + branch
    # t <- t + 1; #remember to comment this
    # print(r); #remember to comment this
  } #lik_nom

  lik <- lik_nom * lik_den

  return(lik)
}

#' Custom likelihood function that applies the procedure after tree splitting. It uses lik_custom for each lineages.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom_split <- function(dataset, N0 = 1, input_check = 1){

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus

  split_matrices <- split_times_matrix(times_matrix)
  Ntips <- length(split_matrices)
  lik   <- rep(NA, Ntips)
  for (i in 1:Ntips)
  {
    ti_i <- split_matrices[[i]][1,1]
    tf_i <- split_matrices[[i]][1, ncol(split_matrices[[i]])]
    tb_i <- NULL
    ts_i <- cbind(split_matrices[[i]][1:2, -c(1, ncol(split_matrices[[i]]))])
    ts_i[2,] <- 1
    if (length(ts_i)==0){ts_i=NULL}
    lambdas_i <- lambdas[split_matrices[[i]][3,]]
    mus_i     <- mus[split_matrices[[i]][3,]]
    tm_i      <- arrange_times_matrix(ti = ti_i, tf = tf_i, tb = tb_i, ts = ts_i)
    dataset_i <- list(times_matrix = tm_i, lambdas = lambdas_i, mus = mus_i)
    lik[i]    <- lik_custom(dataset = dataset_i, N0 = 1)
  }
  return(prod(lik))
}

#' Custom likelihood function for single lineages. You can use it after you apply split_times_matrix.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom_single_lineage   <- function(dataset, N0 = 1, input_check = 1){

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus
  Ntimepoints  <- ncol(times_matrix)
  tf           <- times_matrix2t_coordinates(times_matrix = times_matrix)$tf

  lik_den <- 1; regime <- 1
  if (Ntimepoints > 2){
    for (t in 2:(Ntimepoints - 1))
    {
      time_interval1 <- times_matrix[1,t] - times_matrix[1,t-1]; time_interval2 <- tf - times_matrix[1,t]
      den_term <- lik_single_event(time1 = time_interval1, time2 = time_interval2, lambda = lambdas[regime], mu = mus[regime])
      lik_den  <- lik_den * den_term
      regime   <- regime + 1
    }
  }

  lik_num <- 1; regime <- 1
  for (t in 2:Ntimepoints)
  {
    time_interval <- times_matrix[1,t] - times_matrix[1,t-1]
    num_term <- lik_Pi(lambda = lambdas[regime], mu = mus[regime], time = time_interval)
    lik_num  <- lik_num * num_term
    regime   <- regime + 1
  }

  lik <- unname(lik_num * lik_den)
  return(lik)
}

#' Custom likelihood function that applies the procedure after tree splitting. It uses lik_custom_single_lineage for each lineages.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom_split2 <- function(dataset, N0 = 1, input_check = 1){

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus

  split_matrices <- split_times_matrix(times_matrix)
  Ntips <- length(split_matrices)
  lik   <- rep(NA, Ntips)
  ds <- vector("list", Ntips)
  for (i in 1:Ntips)
  {
    ti_i <- split_matrices[[i]][1,1]
    tf_i <- split_matrices[[i]][1, ncol(split_matrices[[i]])]
    tb_i <- NULL
    ts_i <- cbind(split_matrices[[i]][1:2, -c(1, ncol(split_matrices[[i]]))])
    ts_i[2,] <- 1
    if (length(ts_i)==0){ts_i=NULL}
    lambdas_i <- lambdas[split_matrices[[i]][3,]]
    mus_i     <- mus[split_matrices[[i]][3,]]
    tm_i      <- arrange_times_matrix(ti = ti_i, tf = tf_i, tb = tb_i, ts = ts_i)
    dataset_i <- list(times_matrix = tm_i, lambdas = lambdas_i, mus = mus_i)
    lik[i]    <- lik_custom_single_lineage(dataset = dataset_i, N0 = 1)
  }
  return(prod(lik))
}

#' Custom likelihood function that applies the procedure after tree splitting. It uses lik_custom_single_lineage for each lineages.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom_split3 <- function(dataset, N0 = 1, input_check = 1){
  #uses split_times_matrix3 and lik_custom_single_lineage3

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus

  split_matrices <- split_times_matrix3(times_matrix)
  Ntips <- length(split_matrices)
  lik   <- rep(NA, Ntips)
  for (i in 1:Ntips)
  {
    ti_i <- split_matrices[[i]][1,1]
    tf_i <- split_matrices[[i]][1, ncol(split_matrices[[i]])]
    tb_i <- NULL
    ts_i <- cbind(split_matrices[[i]][1:2, -c(1, ncol(split_matrices[[i]]))])
    ts_i[2,] <- 1
    if (length(ts_i)==0){ts_i=NULL}
    lambdas_i <- lambdas[split_matrices[[i]][3,]]
    mus_i     <- mus[split_matrices[[i]][3,]]
    tm_i      <- arrange_times_matrix(ti = ti_i, tf = tf_i, tb = tb_i, ts = ts_i)
    dataset_i <- list(times_matrix = tm_i, lambdas = lambdas_i, mus = mus_i)
    lik[i]    <- lik_custom_single_lineage3(dataset = dataset_i, N0 = 1)
  }
  return(prod(lik))
}

#' Custom likelihood function for single lineages. You can use it after you apply split_times_matrix.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom_single_lineage3   <- function(dataset, N0 = 1, input_check = 1){

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus
  Ntimepoints  <- ncol(times_matrix)
  tf           <- times_matrix2t_coordinates(times_matrix = times_matrix)$tf

  lik_den <- 1; regime <- 1; shift <- branch <- 0
  if (Ntimepoints > 2){
    for (t in 2:(Ntimepoints - 1))
    {
      shift  <- (sign(times_matrix[2,t]) < 0)
      branch <- (sign(times_matrix[2,t]) > 0)

      time_interval1 <- times_matrix[1,t] - times_matrix[1,t-1]
      time_interval2 <- tf - times_matrix[1,t]
      den_term <- lik_single_event(time1 = time_interval1, time2 = time_interval2, lambda = lambdas[regime], mu = mus[regime])^(1 + branch)
      lik_den  <- lik_den * den_term
      regime   <- regime + shift
    }
  }

  lik_num <- 1; regime <- 1; shift <- branch <- 0
  for (t in 2:Ntimepoints)
  {
    shift  <- (sign(times_matrix[2,t]) < 0)
    branch <- (sign(times_matrix[2,t]) > 0)

    time_interval <- times_matrix[1,t] - times_matrix[1,t-1]
    num_term <- lik_Pi(lambda = lambdas[regime], mu = mus[regime], time = time_interval)
    lik_num  <- lik_num * num_term
    regime   <- regime + shift
  }

  lik <- unname(lik_num * lik_den)
  return(lik)
}
