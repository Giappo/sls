#' RJCB: I think this should move to the tests folder
#' RJCB: I re-use the parameter documentation from
#' one dummy function!
#' @inheritParams default_params_doc
#' @return no idea
#' @export
#' @author Giovanni Laudanno
test_likelihood_formula  <- function(lambdas, mus, ti, tb, ts, tf, N0 = 1, Nsims = 100000,
                           lik_function = lik_custom, sim_function = sim_custom, input_check = 1){

    #the aim is to get "lik_result" equal to "sim_result"
    if (input_check == 1)
    {
      coherent_input <- check_input_data_coherence(lambdas = lambdas, mus = mus, ti = ti, tf = tf, tb = tb, ts = ts, N0 = N0)
      if (coherent_input == 0){stop("Input data are incoherent")}
    }
    time1 <- Sys.time()
    res <- list();  total <- 0; ok <- rep(NA, Nsims)
    while (total < Nsims)
    {
      res    <- sim_function(lambdas = lambdas, mus = mus, ti = ti, tb = tb, ts = ts, tf = tf, input_check = 0); res
      total  <- total + 1
      ok[total] <- res$ok
    }

    #results
    lik_result   <- lik_function(lambdas = lambdas, mus = mus, ti = ti, tb = tb, ts = ts, tf = tf, input_check = 0)
    sim_result   <- sum(ok)/total
    sim_std <- spread <- 0; if (Nsims >= 10000 && sum(ok) >= 10){
      std_results       <- get_std2(oks = ok, lik_result = lik_result, sim_result = sim_result)
      sim_std           <- std_results$std_max
      figure.error_bars <- std_results$figure.error_bars
      spread  <- abs(lik_result - sim_result)/sim_std
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

#' @export
test_likelihood_formula2 <- function(Nsims, dataset, lik_function = lik_custom, sim_function = sim_custom){
  lambdas <- unlist(dataset$lambdas)
  mus     <- unlist(dataset$mus)
  ti      <- unlist(dataset$ti)
  tb      <- unlist(dataset$tb)
  ts      <- unlist(dataset$ts)
  tf      <- unlist(dataset$tf)
  result  <- test_likelihood_formula(lambdas = lambdas, mus = mus,
                                     ti = ti, tb = tb, ts = ts, tf = tf,
                                     Nsims = Nsims, N0 = 1,
                                     lik_function = lik_function, sim_function = sim_function, input_check = 1)
  return(result)
}

# file.edit(paste0(getwd(),"/scripts/main.R"))
