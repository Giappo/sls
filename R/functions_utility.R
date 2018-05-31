# library("testit"); library("DDD"); library("xlsx"); library(MASS)
# general utility functions

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
check_input_data_coherence <- function(lambdas, mus, ti, tf, tb, ts, N0){

  coherent <- 1
  if (length(ts) == 0) {number_of_shifts <- 0} else {number_of_shifts <- ncol(ts)}
  if ( !all.equal(length(lambdas), length(mus), (1+number_of_shifts)) )
  {
    # RJCB: I think you want to use 'stop' or 'warning' here
    stop("Parameters and/or shift times are not coherent.")
  }

  if (is.matrix(tb)){TB <- tb[1,]}else{TB <- tb[1]}
  ntt <- matrix( rbind(c(ti, TB), seq(from = N0, to = length(c(ti,TB)), by = 1 )), nrow=2);
  rownames(ntt) <- c("times", "number")

  if(!is.null(ncol(ts)))
  {
    col_ts  <- list()
    col_ntt <- list()
    for (i in 1:ncol(ts))
    {
      col_ts[[i]] <- ts[,i]
      for (j in 1:ncol(ntt))
      {
        col_ntt[[j]] <- ntt[,j]
        if (ts[1,i] < ntt[1,j]){ coherent <- coherent * (ts[2,i] <= ntt[2,j]) }
      }
    }
  }
  if (coherent == 0)
  {
    cat("Shifting ids are not coherent with branching time")
  }
  return(coherent)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
sim_check_id_presence      <- function(lineages, id){
  ok <- 0
  if (!is.matrix(lineages))
  {
    if (id == lineages[1]){ok <- 1}
  }else
  {
    if (id %in% lineages[,1]){ok <- 1}
  }
  if (id == 0){ok <- 1} #this happens only at the end
  return(ok)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
arrange_times_matrix       <- function(ti, tb, ts, tf){
  ts2  <- ts; if (!is.null(ncol(ts))){ts2 <- -abs(ts)}
  times_matrix  <- cbind(c(ti,0), cbind(ts2, tb), c(tf,0)); times_matrix <- times_matrix[,order(abs(times_matrix[1,]))]; times_matrix[1,] <- abs(times_matrix[1,])
  nshifts <- 0; if(!is.null(ncol(ts))){nshifts <- ncol(ts)}; testit::assert(sum(times_matrix < 0) == nshifts)
  rownames(times_matrix) <- c("when", "who")

  return(times_matrix)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
shuffle                    <- function(oks, N){
  Nmax <- floor(log10(length(oks)))
  Nsamples <- 10^(Nmax - N)
  shuffled_oks <- DDD::sample2(x = oks)
  shuffled_oks <- shuffled_oks[1:(10^Nmax)]
  chunk_size <- 10^N
  subsets <- list(Nsamples)
  subsets <- split(shuffled_oks, ceiling(seq_along(shuffled_oks)/chunk_size))
  medie <- unlist(lapply(X = subsets, FUN = mean))
  return(medie)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
get_std                    <- function(oks){
  oks1 <- oks
  expos <- 2:(floor(log10(length(oks1))) - 1)
  res <- vector("list", (floor((Nmax <- log10(length(oks1)))) - 1))
  fit <- para <- vector("list", (floor(Nmax) - 1))
  for (i in expos){
    res[[i]]  <- shuffle(oks = oks1, N = i)
    fit[[i]]  <- MASS::fitdistr(res[[i]], "normal")
    para[[i]] <- fit[[i]]$estimate
  }
  means <- sds <- rep(NA, length(expos))
  for (i in expos){
    sds[i]   <- para[[i]][2]
    means[i] <- para[[i]][1]
  }
  d.f <- data.frame(
    N    <- expos,
    stds <- log10( sds[-which(is.na(sds))] )
  ); colnames(d.f) <- c("N","std")
  std_fit = stats::lm(stds ~ N, data = d.f)
  x_values <- toString(10^d.f$N)

  ggplot2::ggplot(data = d.f, ggplot2::aes(x = N, y = stds)) + ggplot2::geom_point() + ggplot2::geom_smooth(method='lm') +
    ggplot2::labs(x = "Sample size (log10)", y = "Standard Deviation (log10)")

  log_std_max <- std_fit$coefficients[1] + std_fit$coefficients[2]*Nmax
  std_max <- unname(10^(log_std_max))

  return(list(std_max = std_max, stds = sds, means = means))
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
get_std2                   <- function(oks, lik_result, sim_result){
  repetitions <- 5
  sim_errors <- rep(NA, repetitions)
  stds_vs_n <- means_vs_n <- vector("list", repetitions)
  for (jj in 1:repetitions){
    temp <- sls::get_std(oks)
    sim_errors[jj]   <- temp$std_max
    means_vs_n[[jj]] <- temp$means
    stds_vs_n[[jj]]  <- temp$stds
  }
  sim_std   <- mean(sim_errors)
  mean_vs_n <- colMeans(matrix(unlist(means_vs_n), nrow = repetitions,byrow = T))
  std_vs_n  <- colMeans(matrix(unlist(stds_vs_n ), nrow = repetitions,byrow = T))
  xx <- c(2:(length(mean_vs_n)), log10(length(oks)))
  mean_vs_n <- c(mean_vs_n, sim_result)
  std_vs_n  <- c(std_vs_n , sim_std   )

  sd <- NULL; rm(sd) # nolint, fixes warning: no visible binding for global variable
  para2 <- data.frame( mean = mean_vs_n[-1], sd = std_vs_n[-1] )

  figure <- ggplot2::ggplot(para2, ggplot2::aes(y = mean, x = xx )) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd, ymax = mean + sd)) +
    ggplot2::geom_point() + ggplot2::geom_hline(yintercept = lik_result) +
    ggplot2::xlab("Sample size (log10)") + ggplot2::ylab("Likelihood")

  print(figure)

  return(list(std_max = sim_std, figure.error_bars = figure))
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
check_number_of_species    <- function(lineages){
  ok <- 0
  if (!is.matrix(lineages))
  {
    id <- lineages[1]
    N  <- lineages[2]
    if (id > 0 && N > 0){ok <- 1}
  }else
  {
    ids <- lineages[,1]
    Ns  <- lineages[,2]
    ok  <- prod(N[ids > 0] > 0)
  }
  return(ok)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
export_results_to_xls      <- function(Nsims, sim_function, lik_function, result, sample_threshold = 100000){
  if (Nsims >= sample_threshold && all.equal(sim_function, sim_custom) && all.equal(lik_function, lik_custom))
  {
    results_file      <- paste0(getwd(),"//results//table2.xls")
    figure_name       <- paste0(getwd(),"//results//figura_errori.png")
    figure.error_bars <- result$figure.error_bars
    sheet             <- xlsx::createSheet(results_file, sheetName = result$sheet_name)
    xlsx::write.xlsx(x = result$results.table, file = results_file, sheetName = sheet, append = TRUE)
    # xlsx::addPicture(file = results_file, sheet = test_result$sheet_name, startRow = 13, startColumn = 1)

    # create a png plot
    grDevices::png(figure_name, height=800, width=800, res=250, pointsize=8)
    figure.error_bars
    grDevices::dev.off()
    # Create a new sheet to contain the plot
    # sheet <-createSheet(wb, sheetName = "boxplot")
    # Add title
    # xlsx.addTitle(sheet, rowIndex=1, title="Error as function of sample size",
    #               titleStyle = TITLE_STYLE)
    # Add the plot created previously
    xlsx::addPicture(file = figure_name, sheet = sheet, scale = 1, startRow = 13, startColumn = 1)
    # remove the plot from the disk
    res <- file.remove(figure_name)
    # Save the workbook to a file...
    #++++++++++++++++++++++++++++++++++++
    # saveWorkbook(wb, "r-xlsx-report-example.xlsx")
  }
}
#install.packages(r2excel)

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
export_results_to_xls_2   <- function(Nsims, sim_function, lik_function, result, sample_threshold = 100000){
  if (Nsims >= sample_threshold && all.equal(sim_function, sim_custom) && all.equal(lik_function, lik_custom))
  {
    results_file      <- paste0(getwd(),"//results//table2.xlsx")
    figure_name       <- paste0(getwd(),"//results//figura_errori.png")
    figure.error_bars <- result$figure.error_bars

    xlsx::write.xlsx(x = result$results.table, file = results_file, sheetName = result$sheet_name, append = TRUE)

    wb    <- xlsx::createWorkbook()
    sheet <- xlsx::createSheet(wb, result$sheet_name)
    d.f   <- as.data.frame(result$results.table)

    grDevices::png(figure_name, height=1000, width=1000, res=300, pointsize=8)
    figure.error_bars
    grDevices::dev.off()

    startRow = nrow(result$results.table) + 2
    xlsx::addPicture(file = figure_name, sheet = sheet, scale = 1, startRow = startRow,
               startColumn = 1)


    # xlsx::addDataFrame(d.f, sheet = sheet, row.names=FALSE, col.names=FALSE, startRow = 1, showNA = F)
    xlsx::saveWorkbook(wb, results_file)


  }
}

#' Determines if the environment is Travis CI
#' @return TRUE if run on Travis CI, FALSE otherwise
#' @author Richel J.C. Bilderbeek
#' @export
is_on_travis <- function() {
  Sys.getenv("TRAVIS") != ""
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
split_times_matrix <- function(times_matrix, N0 = 1){

  # ###
  # ti <- d.s$ti; tb <- d.s$tb; ts <- d.s$ts; tf <- d.s$tf; N0 <- 1
  # times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
  # ###

  #determine rates
  # tb <- cbind(times_matrix[1:2, which(sign(times_matrix[,2])>0)])
  tb <- times_matrix2t_coordinates(times_matrix = times_matrix)$tb
  ti <- times_matrix2t_coordinates(times_matrix = times_matrix)$ti
  nbranches <- 0; if(!is.null(ncol(tb))){nbranches <- ncol(tb)};
  Ntips <- N0 + nbranches; r <- rep(regime <- 1, Ntips); n <- N0
  rr <- vector("list", ncol(times_matrix)); rr[[1]] <- r
  for (t in 2:ncol(times_matrix))
  {
    upds <- update_regimes(t = t, regime = regime, times_matrix = times_matrix, r = r, n = n)
    rr[[t]] <- r <- upds$r
    regime  <- upds$regime
    n = n + (sign(times_matrix[2,t])==1)
  }
  rrr <- matrix(unlist(rr), nrow = Ntips)
  times_matrix1 <- rbind(times_matrix, rrr);

  #TM is the list of splitted matrices
  TM <- vector("list",Ntips)
  birth_times <- c(ti, tb[1,])
  all_times <- times_matrix1[1,]
  tb_positions <- match(birth_times, all_times)
  for (i in 1:Ntips)
  {
    TM[[i]] <- times_matrix1[c(1,2,2+i),c(tb_positions[i], which(times_matrix1[2,]==-i), ncol(times_matrix1))]
    TM[[i]][2,1] <- 0
    TM[[i]][3,ncol(TM[[i]])] <- 0
    rownames(TM[[i]])[3] <- "regime"
  };TM
  return(TM)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
split_times_matrix3 <- function(times_matrix, N0 = 1){

  # ###
  # ti <- d.s$ti; tb <- d.s$tb; ts <- d.s$ts; tf <- d.s$tf; N0 <- 1
  # times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
  # ###

  #determine rates
  # tb <- cbind(times_matrix[1:2, which(sign(times_matrix[,2])>0)])
  tb <- times_matrix2t_coordinates(times_matrix = times_matrix)$tb
  ti <- times_matrix2t_coordinates(times_matrix = times_matrix)$ti
  nbranches <- 0; if(!is.null(ncol(tb))){nbranches <- ncol(tb)};
  Ntips <- N0 + nbranches; r <- rep(regime <- 1, Ntips); n <- N0
  rr <- vector("list", ncol(times_matrix)); rr[[1]] <- r
  for (t in 2:ncol(times_matrix))
  {
    upds <- update_regimes(t = t, regime = regime, times_matrix = times_matrix, r = r, n = n)
    rr[[t]] <- r <- upds$r
    regime  <- upds$regime
    n = n + (sign(times_matrix[2,t])==1)
  }
  rrr <- matrix(unlist(rr), nrow = Ntips)
  times_matrix1 <- rbind(times_matrix, rrr);

  #TM is the list of splitted matrices
  TM <- vector("list",Ntips)
  birth_times <- c(ti, tb[1,])
  all_times <- times_matrix1[1,]
  tb_positions <- match(birth_times, all_times)
  condx <- function(times_matrix1, i, tb_positions)
  {
    x <- times_matrix1[1,]
    y <- times_matrix1[2,]
    (x >= x[tb_positions[i]]) & ((y == -i) | (y > 0) | y == y[length(y)])
  }
  for (i in 1:Ntips)
  {
    TM[[i]] <- times_matrix1[c(1,2,2+i), condx(times_matrix1 = times_matrix1, i = i, tb_positions = tb_positions)]
    TM[[i]][2,1] <- 0
    TM[[i]][3,ncol(TM[[i]])] <- 0
    rownames(TM[[i]])[3] <- "regime"
  };TM
  return(TM)
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
update_regimes <- function(t, regime, times_matrix, r, n, N0 = 1){

  # nbranches <- 0; if(!is.null(ncol(tb))){nbranches <- ncol(tb)};

  Ntips <- N0 + sum(sign(times_matrix[2,])>0)

  if (sign(times_matrix[2,t])!=-1){return(list(r = r, regime = regime))}
  ttb <- times_matrix[2,sign(times_matrix[2,]) == +1]
  fathers <- unname( c(rep(0, N0), ttb) ); sons <- 1:Ntips
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
  return(list(r = r, regime = regime))
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
times_matrix2t_coordinates <- function(times_matrix){
  ti <- unname(times_matrix[1,1])
  tf <- unname(times_matrix[1, ncol(times_matrix)])
  tb <- cbind(times_matrix[,which(times_matrix[2,] > 0)])     ; if(!is.null(tb)){rownames(tb) <- c("tb_when","tb_who")}
  ts <- abs(cbind(times_matrix[,which(times_matrix[2,] < 0)])); if(!is.null(ts)){rownames(ts) <- c("ts_when","ts_who")}
  if (sum(ts) == 0 && prod(ts) == 1){ts = NULL}
  if (sum(tb) == 0 && prod(tb) == 1){tb = NULL}
  # ts_i <- cbind((times_matrix[1:2, -c(1,ncol(times_matrix))]))
  return(list(ti = ti, tf = tf, ts = ts, tb = tb))
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
variable_name2string <- function(v1){
  deparse(substitute(v1))
}

#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
load_all_data <- function(the.environment = environment()){
  d <- data(package = "sls",envir = the.environment)
  d$results[, "Item"]
  nm <- d$results[, "Item"]
  data(list = nm, package = "sls",envir = the.environment)
}




