parts_lik_custom      <- function(dataset, N0 = 1, input_check = 1){

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
  Ntips   <-  N0 + nbranches
  fathers <- unname( c(rep(0, N0), tb[2,]) ); sons <- 1:Ntips
  lik_den <- 1; r <- rep(regime <- 1, Ntips); n <- 1; t <- 2; shift <- branch <- 0
  parts_den <- list(); p <- 1
  for (t in 2:(Ntimepoints - 1))
  {
    shift  <- (sign(times_matrix[2,t]) < 0)
    branch <- (sign(times_matrix[2,t]) > 0)

    lineages_affected_by_the_event <- rep(0,n); lineages_affected_by_the_event[abs(times_matrix[2,t])] <- 1
    expo <- unname(2 - lineages_affected_by_the_event * shift)
    # print(lineages_affected_by_the_event)

    time_interval1 <- times_matrix[1,t] - times_matrix[1,t-1]; time_interval2 <- tf - times_matrix[1,t]
    per_lineage_outcome <- lik_single_event(time1 = time_interval1, time2 = time_interval2, lambda = lambdas[r[1:n]], mu = mus[r[1:n]])^expo
    parts_den[[p]] <- per_lineage_outcome; p <- p + 1
    lik_den <- lik_den * prod(per_lineage_outcome)

    if (shift)
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
  parts_num <- list(); p <- 1
  for (t in 2:Ntimepoints)
  {
    shift  <- (sign(times_matrix[2,t]) < 0)
    branch <- (sign(times_matrix[2,t]) > 0)

    # print(r[1:n])
    # print(lambdas[r[1:n]])
    time_interval <- times_matrix[1,t] - times_matrix[1,t-1]
    per_lineage_outcome <- lik_Pi(lambda = lambdas[r[1:n]], mu = mus[r[1:n]], time = time_interval)
    parts_num[[p]] <- per_lineage_outcome; p <- p + 1
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

  return(list(lik = lik, parts_num = parts_num, parts_den = parts_den))
}

parts_lik_custom_split2 <- function(dataset, N0 = 1, input_check = 1){

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus

  split_matrices <- split_times_matrix(times_matrix)
  Ntips <- length(split_matrices)
  # lik   <- rep(NA, Ntips)
  parts <- parts_num <- parts_den <- list()
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
    parts     <- parts_lik_custom_single_lineage(dataset = dataset_i, N0 = 1)
    parts_num[[i]] <- unlist(parts$parts_num)
    parts_den[[i]] <- unlist(parts$parts_den)
  }
  return(list(parts_num = parts_num, parts_den = parts_den))
}

parts_lik_custom_single_lineage   <- function(dataset, N0 = 1, input_check = 1){

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus
  Ntimepoints  <- ncol(times_matrix)
  tf           <- times_matrix2t_coordinates(times_matrix = times_matrix)$tf

  lik_den <- 1; regime <- 1
  parts_den <- list(); p <- 1
  if (Ntimepoints > 2){
    for (t in 2:(Ntimepoints - 1))
    {
      time_interval1 <- times_matrix[1,t] - times_matrix[1,t-1]; time_interval2 <- tf - times_matrix[1,t]
      den_term <- lik_single_event(time1 = time_interval1, time2 = time_interval2, lambda = lambdas[regime], mu = mus[regime])
      parts_den[[p]] <- den_term; p <- p + 1
      lik_den  <- lik_den * den_term
      regime   <- regime + 1
    }
  }else{parts_den <- 1}

  lik_num <- 1; regime <- 1
  parts_num <- list(); p <- 1
  for (t in 2:Ntimepoints)
  {
    time_interval <- times_matrix[1,t] - times_matrix[1,t-1]
    num_term <- lik_Pi(lambda = lambdas[regime], mu = mus[regime], time = time_interval)
    parts_num[[p]] <- num_term; p <- p + 1
    lik_num  <- lik_num * num_term
    regime   <- regime + 1
  }

  lik <- unname(lik_num * lik_den)
  return(list(lik = lik, parts_num = parts_num, parts_den = parts_den))
}

####
dataset <- dataset_Bart

test1 <- list()
temp <- parts_lik_custom(dataset)
test1$num <- temp$parts_num
test1$den <- temp$parts_den

test2 <- list()
temp <- parts_lik_custom_split2(dataset)
test2$num <- temp$parts_num
test2$den <- temp$parts_den

print(test1)
print(test2)
dataset

prod(unlist(test1$num)) * prod(unlist(test1$den))
prod(unlist(test2$num)) * prod(unlist(test2$den))


test2_branch1 <- prod(test2$num[[1]]) * prod(test2$den[[1]])
test2_branch2 <- prod(test2$num[[2]]) * prod(test2$den[[2]])

branch1_test1_num <- 1; for (i in 1: length(test1$num)) {branch1_test1_num <- branch1_test1_num * test1$num[[i]][1]}
branch1_test1_den <- 1; for (i in 1: length(test1$den)) {branch1_test1_den <- branch1_test1_den * test1$den[[i]][1]}

branch2_test1_num <- 1
for (i in 1: length(test1$num)){
  temp <- test1$num[[i]][2]; if(is.na(temp)){temp <- 1}
  branch2_test1_num <- branch2_test1_num * temp
}
branch2_test1_den <- 1
for (i in 1: length(test1$den)){
  temp <- test1$den[[i]][2]; if(is.na(temp)){temp <- 1}
  branch2_test1_den <- branch2_test1_den * temp
}

branch1_test1_num
branch1_test1_den

branch2_test1_num
branch2_test1_den

test1_branch1 = branch1_test1_num * branch1_test1_den
test1_branch2 = branch2_test1_num * branch2_test1_den

test1_branch1
test1_branch2
test2_branch1
test2_branch2


t1_branch1_num_list <- list(); for (i in 1: length(test1$num)) {t1_branch1_num_list[[i]] <- test1$num[[i]][1]}
t1_branch1_den_list <- list(); for (i in 1: length(test1$den)) {t1_branch1_den_list[[i]] <- test1$den[[i]][1]}

t2_branch1_num_list <- list(); for (i in 1: length(test2$num[[1]])) {t2_branch1_num_list[[i]] <- test2$num[[1]][i]}
t2_branch1_den_list <- list(); for (i in 1: length(test2$den[[1]])) {t2_branch1_den_list[[i]] <- test2$den[[1]][i]}

t1_b1_num <- unname(unlist(t1_branch1_num_list));t1_b1_num
t1_b1_den <- unname(unlist(t1_branch1_den_list));t1_b1_den
t2_b1_num <- unname(unlist(t2_branch1_num_list));t2_b1_num
t2_b1_den <- unname(unlist(t2_branch1_den_list));t2_b1_den

fragments <- function (t1_b1_num, t1_b1_den){
t1_b1 <- rep(NA, length(t1_b1_den)); for(i in 1:length(t1_b1_den)){t1_b1[i] <- t1_b1_num[i] * t1_b1_num[i+1] * t1_b1_den[i]}
return(t1_b1)
}
t1_frags <- fragments(t1_b1_num, t1_b1_den)
t2_frags <- fragments(t2_b1_num, t2_b1_den)

t1_frags[1] * t1_frags[2]
t1_frags[2] * t1_frags[3]
t2_frags[1] * t2_frags[2]
