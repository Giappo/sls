#' @inheritParams default_params_doc
#' @return no idea
#' @export
#' @author Giovanni Laudanno
test_likelihood_formula  <- function(dataset, N0 = 1, Nsims = 100000,
                           lik_function = lik_custom, sim_function = sim_custom, input_check = TRUE){

  times_matrix <- dataset$times_matrix
  lambdas      <- dataset$lambdas
  mus          <- dataset$mus
  coords       <- times_matrix2t_coordinates(times_matrix = times_matrix)
  ti           <- coords$ti
  tb           <- coords$tb
  ts           <- coords$ts
  tf           <- coords$tf

  #the aim is to get "lik_result" equal to "sim_result"
  if (input_check == 1)
  {
    coherent_input <- check_input_data_coherence(dataset = dataset, N0 = N0)
    if (coherent_input == 0){stop("Input data are incoherent")}
  }
  time1 <- Sys.time()
  # res <- list();  total <- 0; ok <- rep(NA, Nsims)
  # while (total < Nsims)
  # {
  #   res    <- sim_function(dataset = dataset, input_check = 0); res
  #   total  <- total + 1
  #   ok[total] <- res$ok
  # }

  #results
  lik_result   <- lik_function(dataset = dataset, input_check = 0)
  sim_out      <- sim_custom_repeat(dataset = dataset, Nsims = Nsims, sim_function = sim_function, N0 = N0, input_check = FALSE)
  sim_result   <- sim_out$sim_result
  ok           <- sim_out$ok

  figure.error_bars <- NULL; sim_std <- spread <- 0; if (Nsims >= 10000 && sum(ok) >= 10){
    std_results       <- get_std2(oks = ok, lik_result = lik_result, sim_result = sim_result)
    sim_std           <- std_results$std_max
    figure.error_bars <- std_results$figure.error_bars
    spread  <- (lik_result - sim_result)/sim_std
    cat(paste0("sim result and lik result agree within ", signif(x = spread, digits = 2) ," standard deviations.\n"))
  }
  time_elapsed <- as.double(difftime(Sys.time(), time1, units = "secs"))

  #saveable output
  datas <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
  results <- matrix(NA, nrow = nrow(datas) + 4, ncol = ncol(datas))
  results[,1] <- c(lik_result, sim_result, sim_std, spread, Nsims, time_elapsed);

  rates <- matrix(NA, nrow = 2, ncol = ncol(datas))
  rates[1, datas[2,] < 0 | datas[1,] == datas[1, 1]] <- lambdas
  rates[2, datas[2,] < 0 | datas[1,] == datas[1, 1]] <- mus

  results.table <- rbind(datas, rates, results)
  rownames(results.table) <- c("when", "who", "lambdas", "mus", "lik_results", "sim_results", "sim_std", "spread", "Nsims", "time(seconds)")
  sheet_name <- toString(format(Sys.time(), "%Y-%m-%d %H.%M.%S") )

  return(list(lik_res = lik_result, sim_res = sim_result, sim_std = sim_std, spread = spread, results.table = results.table, sheet_name = sheet_name, figure.error_bars = figure.error_bars, ok_vec = ok))
}


#' Does something
#' @inheritParams default_params_doc
#' @return result
#' @export
test_likelihood_formula2 <- function(s,Nsims){
  #lik_function = lik_custom, sim_function = sim_custom
  load_all_data(the.environment = environment()); data.sets <- ls(pattern = "dataset_",envir = environment())
  dataset <- get(data.sets[[s]])
  test_result <- test_likelihood_formula(dataset = dataset, Nsims = Nsims,
                          sim_function = sim_custom2, lik_function = lik_custom)

  if (Nsims >= 1E5)
  {
    results_file <- paste0(getwd(),"//results//table3.xls")
    sheet_name <- sheet_name0 <- data.sets[[s]]
    wb <- xlsx::loadWorkbook(file = results_file)
    prova <- xlsx::getSheets(wb = wb)
    nomi <- names(prova)
    jj <- 1; while (sheet_name %in% nomi) {jj <- jj + 1; sheet_name <- paste0(sheet_name0," - ", toString(jj))}

    xlsx::write.xlsx(x = test_result$results.table, file = results_file, sheetName = sheet_name, append = TRUE)
  }
  return(test_result)
}
