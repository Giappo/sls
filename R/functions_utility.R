# library("testit"); library("DDD"); library("xlsx"); library(MASS)
# general utility functions
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
#' @export
arrange_times_matrix       <- function(ti, tb, ts, tf){
  ts2  <- ts; if (!is.null(ncol(ts))){ts2 <- -ts}
  times_matrix  <- cbind(c(ti,0), cbind(ts2, tb), c(tf,0)); times_matrix <- times_matrix[,order(abs(times_matrix[1,]))]; times_matrix[1,] <- abs(times_matrix[1,])
  nshifts <- 0; if(!is.null(ncol(ts))){nshifts <- ncol(ts)}; testit::assert(sum(times_matrix < 0) == nshifts)
  rownames(times_matrix) <- c("when", "who")

  return(times_matrix)
}
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
#' @export
get_std2                   <- function(oks, lik_result, sim_result){
  repetitions <- 5
  sim_errors <- rep(NA, repetitions)
  stds_vs_n <- means_vs_n <- vector("list", repetitions)
  for (jj in 1:repetitions){
    temp <- sls:::get_std(oks)
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
#' @export
export_results_to_xls      <- function(Nsims, sim_function, lik_function, result, sample_threshold = 100000){
  if (Nsims >= sample_threshold && all.equal(sim_function, sim_custom) && all.equal(lik_function, lik_custom))
  {
    results_file      <- paste0(getwd(),"//results//table2.xls")
    figure_name       <- paste0(getwd(),"//results//figura_errori.png")
    figure.error_bars <- result$figure.error_bars
    sheet             <- xlsx::createSheet(results_file, sheetName = result$sheet_name)
    xlsx:::write.xlsx(x = result$results.table, file = results_file, sheetName = sheet, append = TRUE)
    # xlsx:::addPicture(file = results_file, sheet = test_result$sheet_name, startRow = 13, startColumn = 1)

    # create a png plot
    png(figure_name, height=800, width=800, res=250, pointsize=8)
    figure.error_bars
    dev.off()
    # Create a new sheet to contain the plot
    # sheet <-createSheet(wb, sheetName = "boxplot")
    # Add title
    # xlsx.addTitle(sheet, rowIndex=1, title="Error as function of sample size",
    #               titleStyle = TITLE_STYLE)
    # Add the plot created previously
    xlsx:::addPicture(file = figure_name, sheet = sheet, scale = 1, startRow = 13, startColumn = 1)
    # remove the plot from the disk
    res <- file.remove(figure_name)
    # Save the workbook to a file...
    #++++++++++++++++++++++++++++++++++++
    # saveWorkbook(wb, "r-xlsx-report-example.xlsx")
  }
}
install.packages(r2excel)
export_results_to_xls      <- function(Nsims, sim_function, lik_function, result, sample_threshold = 100000){
  if (Nsims >= sample_threshold && all.equal(sim_function, sim_custom) && all.equal(lik_function, lik_custom))
  {
    results_file      <- paste0(getwd(),"//results//table2.xlsx")
    figure_name       <- paste0(getwd(),"//results//figura_errori.png")
    figure.error_bars <- result$figure.error_bars

    xlsx::write.xlsx(x = result$results.table, file = results_file, sheetName = result$sheet_name, append = TRUE)

    wb    <- xlsx::createWorkbook()
    sheet <- xlsx::createSheet(wb, result$sheet_name)
    d.f   <- as.data.frame(result$results.table)

    png(figure_name, height=1000, width=1000, res=300, pointsize=8)
    figure.error_bars
    dev.off()

    startRow = nrow(result$results.table) + 2
    xlsx::addPicture(file = figure_name, sheet = sheet, scale = 1, startRow = startRow,
               startColumn = 1)


    # xlsx::addDataFrame(d.f, sheet = sheet, row.names=FALSE, col.names=FALSE, startRow = 1, showNA = F)
    xlsx::saveWorkbook(wb, results_file)


  }
}
