#basic functions

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

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_Pix <- function (lambdas, mus, times){
  N  <- length(lambdas)
  Ti <- cumsum((times))
  x  <- rev(Ti[N] - sort(Ti,decreasing = T))

  Pis   <- lik_Pi(lambda = lambdas[1:N], mu = mus[1:N], time = times[1:N])
  us    <- lik_ut(lambda = lambdas[1:(N-1)], mu = mus[1:(N-1)], time = times[1:(N-1)])
  ps    <- lik_pt(lambda = lambdas[1:(N-1)], mu = mus[1:(N-1)], time = x[1:(N-1)])

  likelihood <- prod( Pis ) * ( prod( (1 - us * (1 - ps) ) ) )^-1

  return(likelihood)
}

#likelihoods

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_single_event <- function(lambda, mu, time1, time2){ #if squared describes the linkage between two consecutives full lines, if its exponent is one it describes instead a single-lineage shift.

  lik <- 1 - lik_ut(lambda = lambda, mu = mu, time = time1) *
    (1 - lik_pt(lambda = lambda, mu = mu, time = time2))

  lik <- lik^(-1)

  return(lik)
} #basic denominator likelihood module. if squared yields branching contribution

#' Custom likelihood function.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom      <- function(lambdas, mus, ti, tf, tb, ts, N0 = 1, input_check = 1){
  # lik_custom      <- function(lambdas, mus, times_matrix, N0 = 1, input_check = 1){

  if (input_check == 1){
    coherent_input <- check_input_data_coherence(lambdas = lambdas, mus = mus, ti = ti, tf = tf, tb = tb, ts = ts, N0 = N0)
    if (coherent_input == 0){stop("Input data are incoherent")}
  }

  if (!is.null(ncol(tb))){ rownames(tb) <- c("time", "who") }
  if (!is.null(ncol(ts))){ rownames(ts) <- c("time", "who") }
  nbranches <- 0; if(!is.null(ncol(tb))){nbranches <- ncol(tb)};
  Ntips <-  N0 + nbranches
  times_matrix  <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
  fathers <- unname( c(rep(0, N0), tb[2,]) ); sons <- 1:Ntips
  lik_den <- 1; r <- rep(regime <- 1, Ntips); n <- 1; t <- 2; shift <- branch <- 0
  for (t in 2:(ncol(times_matrix) - 1))
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
  for (t in 2:ncol(times_matrix))
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
lik_custom_split <- function(lambdas, mus, ti, tf, tb, ts, N0 = 1, input_check = 1){
  times_matrix   <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
  split_matrices <- split_times_matrix(times_matrix)
  Ntips <- length(split_matrices)
  lik   <- rep(NA, Ntips)
  for (i in 1:Ntips)
  {
    ti_i <- split_matrices[[i]][1,1]
    tf_i <- split_matrices[[i]][1, ncol(split_matrices[[i]])]
    tb_i <- NULL
    ts_i <- cbind(split_matrices[[i]][1:2, -c(1,ncol(split_matrices[[i]]))])
    ts_i[2,] <- 1
    if (length(ts_i)==0){ts_i=NULL}
    lambdas_i <- lambdas[split_matrices[[i]][3,]]
    mus_i     <- mus[split_matrices[[i]][3,]]
    lik[i]    <- lik_custom(lambdas = lambdas_i, mus = mus_i, ti = ti_i, tf = tf_i, tb = tb_i, ts = ts_i)
  }
  return(prod(lik))
}

#' Custom likelihood function for single lineages. You can use it after you apply split_times_matrix.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom_single_lineage   <- function(lambdas, mus, ti, tf, tb, ts, N0 = 1, input_check = 1){
  times_matrix  <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
  Ntimepoints   <- ncol(times_matrix)

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
lik_custom_split2 <- function(lambdas, mus, ti, tf, tb, ts, N0 = 1, input_check = 1){
  times_matrix   <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
  split_matrices <- split_times_matrix(times_matrix)
  Ntips <- length(split_matrices)
  lik   <- rep(NA, Ntips)
  for (i in 1:Ntips)
  {
    ti_i <- split_matrices[[i]][1,1]
    tf_i <- split_matrices[[i]][1, ncol(split_matrices[[i]])]
    tb_i <- NULL
    ts_i <- cbind(split_matrices[[i]][1:2, -c(1,ncol(split_matrices[[i]]))])
    ts_i[2,] <- 1
    if (length(ts_i)==0){ts_i=NULL}
    lambdas_i <- lambdas[split_matrices[[i]][3,]]
    mus_i     <- mus[split_matrices[[i]][3,]]
    lik[i]    <- lik_custom_single_lineage(lambdas = lambdas_i, mus = mus_i, ti = ti_i, tf = tf_i, tb = tb_i, ts = ts_i)
  }
  return(prod(lik))
}

#' Custom likelihood function that applies the procedure after tree splitting. It uses lik_custom_single_lineage for each lineages.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom_split3 <- function(lambdas, mus, ti, tf, tb, ts, N0 = 1, input_check = 1){
  #uses split_times_matrix3 and lik_custom_single_lineage3
  times_matrix   <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
  split_matrices <- split_times_matrix3(times_matrix)
  Ntips <- length(split_matrices)
  lik   <- rep(NA, Ntips)
  for (i in 1:Ntips)
  {
    ti_i <- split_matrices[[i]][1,1]
    tf_i <- split_matrices[[i]][1, ncol(split_matrices[[i]])]
    tb_i <- NULL
    ts_i <- cbind(split_matrices[[i]][1:2, -c(1,ncol(split_matrices[[i]]))])
    ts_i[2,] <- 1
    if (length(ts_i)==0){ts_i=NULL}
    lambdas_i <- lambdas[split_matrices[[i]][3,]]
    mus_i     <- mus[split_matrices[[i]][3,]]
    lik[i]    <- lik_custom_single_lineage3(lambdas = lambdas_i, mus = mus_i, ti = ti_i, tf = tf_i, tb = tb_i, ts = ts_i)
  }
  return(prod(lik))
}

#' Custom likelihood function for single lineages. You can use it after you apply split_times_matrix.
#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
lik_custom_single_lineage3   <- function(lambdas, mus, ti, tf, tb, ts, N0 = 1, input_check = 1){
  times_matrix  <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
  Ntimepoints   <- ncol(times_matrix)

  lik_den <- 1; regime <- 1; shift <- branch <- 0
  if (Ntimepoints > 2){
    for (t in 2:(Ntimepoints - 1))
    {
      time_interval1 <- times_matrix[1,t] - times_matrix[1,t-1]; time_interval2 <- tf - times_matrix[1,t]
      shift  <- (sign(times_matrix[2,t]) < 0)
      branch <- (sign(times_matrix[2,t]) > 0)

      den_term <- lik_single_event(time1 = time_interval1, time2 = time_interval2, lambda = lambdas[regime], mu = mus[regime])^(1 + branch)
      lik_den  <- lik_den * den_term
      regime   <- regime + shift
    }
  }

  lik_num <- 1; regime <- 1; shift <- branch <- 0
  for (t in 2:Ntimepoints)
  {
    time_interval <- times_matrix[1,t] - times_matrix[1,t-1]
    shift  <- (sign(times_matrix[2,t]) < 0)
    branch <- (sign(times_matrix[2,t]) > 0)

    num_term <- lik_Pi(lambda = lambdas[regime], mu = mus[regime], time = time_interval)
    lik_num  <- lik_num * num_term
    regime   <- regime + shift
  }

  lik <- unname(lik_num * lik_den)
  return(lik)
}
