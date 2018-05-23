#basic functions
#' @export
lik_pt  <- function (lambda, mu, time){
  Lambda <- exp((mu-lambda)*(time))
  out    <- (lambda==mu)*( 1/(1+lambda*(time)) ) +
    (lambda!=mu)*( (lambda-mu)/(lambda-mu*Lambda*(lambda!=mu)) )
  return(out)
}
#' @export
lik_ut  <- function (lambda, mu, time){
  Lambda <- exp((mu-lambda)*(time))
  out    <- (lambda==mu)*( lambda*(time)/(1+lambda*(time)) ) +
    (lambda!=mu)*( (lambda - lambda*Lambda) / (lambda - mu*Lambda*(lambda!=mu)) )
  return(out)
}
#' @export
lik_Pi  <- function (lambda, mu, time){
  out <- lik_pt(lambda = lambda, mu = mu, time = time)*(1-lik_ut(lambda = lambda, mu = mu, time = time))
  return(out)
}
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
#' @export
lik_single_event <- function(lambda, mu, time1, time2){ #if squared describes the linkage between two consecutives full lines, if its exponent is one it describes instead a single-lineage shift.

  lik <- 1 - lik_ut(lambda = lambda, mu = mu, time = time1) *
    (1 - lik_pt(lambda = lambda, mu = mu, time = time2))

  lik <- lik^(-1)

  return(lik)
} #basic denominator likelihood module. if squared yields branching contribution

#' Custom likelihood equation (for lack of better words)
#' @export
lik_custom       <- function(lambdas, mus, ti, tf, tb, ts, N0 = 1, input_check = 1){

  if (input_check == 1){
    coherent_input <- check_input_data_coherence(ti = ti, tf = tf, tb = tb, ts = ts)
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
      regime <- regime + 1
      who_shifts <- abs(times_matrix[2,t])
      r[who_shifts] <- regime
      future <- (Ntips > (n)) * ((n+1):Ntips)

      who_shifts_cascade3 <- NULL
      who_shifts_cascade  <- who_shifts
      who_shifts_cascade2 <- sons[which(fathers %in% who_shifts_cascade)]
      who_shifts_cascade3 <- unique(c(who_shifts_cascade, who_shifts_cascade2))
      who_shifts_cascade3 <- intersect(who_shifts_cascade3, future); who_shifts_cascade3
      while ( prod(all.equal(who_shifts_cascade, who_shifts_cascade3) != 1) )
      {
        who_shifts_cascade  <- who_shifts_cascade3
        who_shifts_cascade2 <- sons[which(fathers %in% who_shifts_cascade)]
        who_shifts_cascade3 <- unique(c(who_shifts_cascade, who_shifts_cascade2));
        who_shifts_cascade3 <- intersect(who_shifts_cascade3, future); #print(who_shifts_cascade3)
      }
      who_shifts_cascade    <- who_shifts_cascade3; who_shifts_cascade
      r[who_shifts_cascade] <- r[who_shifts]
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
      regime <- regime + 1
      who_shifts <- abs(times_matrix[2,t])
      r[who_shifts] <- regime
      future <- (Ntips > (n)) * ((n+1):Ntips)

      who_shifts_cascade3 <- NULL
      who_shifts_cascade  <- who_shifts
      who_shifts_cascade2 <- sons[which(fathers %in% who_shifts_cascade)]
      who_shifts_cascade3 <- unique(c(who_shifts_cascade, who_shifts_cascade2))
      who_shifts_cascade3 <- intersect(who_shifts_cascade3, future); who_shifts_cascade3
      while ( prod(all.equal(who_shifts_cascade, who_shifts_cascade3) != 1) )
      {
        who_shifts_cascade  <- who_shifts_cascade3
        who_shifts_cascade2 <- sons[which(fathers %in% who_shifts_cascade)]
        who_shifts_cascade3 <- unique(c(who_shifts_cascade, who_shifts_cascade2));
        who_shifts_cascade3 <- intersect(who_shifts_cascade3, future);# print(who_shifts_cascade3)
      }
      who_shifts_cascade    <- who_shifts_cascade3; who_shifts_cascade
      r[who_shifts_cascade] <- r[who_shifts]
    }
    n <- n + branch
    # t <- t + 1; #remember to comment this
    # print(r); #remember to comment this
  } #lik_nom

  lik <- lik_nom * lik_den

  return(lik)
} #custom likelihood function: requires topology information. requires lik_event_shift
