rm(list = ls()); library(sls)
load_all_data(the.environment = environment()); data.sets <- ls(pattern = "dataset_",envir = environment())
sim_function = sim_custom; lik_function = lik_custom_split2
Nsims <- 1E3 #1E6
start_coord <- 11 #use it in case you already have some partial results
for (i in start_coord:length(data.sets))
{
  print(i)
  d.s <- get(data.sets[[i]])
  #the aim is to get "lik_result" equal to "sim_result" for an high enough number of simulations
  test_result1 <- test_likelihood_formula(dataset = d.s, Nsims = Nsims, lik_function = lik_function, sim_function = sim_function); print(test_result1$results.table)
  if (Nsims >= 100000)
  {
    results_file <- paste0(getwd(),"//results//table2.xls")
    xlsx::write.xlsx(x = test_result1$results.table, file = results_file, sheetName = data.sets[[i]], append = TRUE)
    # xlsx::addPicture(file = results_file, sheet = test_result$sheet_name, startRow = 13, startColumn = 1)
  }
}




# #####-test splits-#####
# rm(list = ls())
# library(sls)
# load_all_data(the.environment = environment())
# data.sets <- ls(pattern = "dataset_", envir = environment())
#
# out1 <- out2 <- out3 <- rep(NA, length(data.sets)); i <- 1
# for (i in 1:length(data.sets))
# {
#   data    <- get(data.sets[[i]])
#   lambdas <- data$lambdas
#   mus     <- data$mus
#   ti      <- data$ti
#   tb      <- data$tb
#   ts      <- data$ts
#   tf      <- data$tf
#   out1[i] <- lik_custom_split(lambdas = lambdas, mus = mus, ti = ti, tb = tb, ts = ts, tf = tf)
#   out2[i] <- lik_custom_split2(lambdas = lambdas, mus = mus, ti = ti, tb = tb, ts = ts, tf = tf)
#   out3[i] <- lik_custom_split3(lambdas = lambdas, mus = mus, ti = ti, tb = tb, ts = ts, tf = tf)
# }
# print(out1)
# print(out2)
# print(out3)
#

####
# rm(list = ls())
# library(sls)
# load_all_data(the.environment = environment())
# # data.sets <- ls(pattern = "dataset_",envir = environment())
# data <- get("dataset_pure_shifting3")
# lambdas <- data$lambdas
# mus <- data$mus
# ti <- data$ti
# tb <- data$tb
# ts <- data$ts
# tf <- data$tf
#
# out1 <- lik_custom(lambdas = lambdas, mus = mus, ti = ti, tb = tb, ts = ts, tf = tf)
# out2 <- lik_custom_single_lineage(lambdas = lambdas, mus = mus, ti = ti, tb = tb, ts = ts, tf = tf)
# print(out1)
# print(out2)














