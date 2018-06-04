rm(list = ls()); library(sls)
load_all_data(the.environment = environment()); data.sets <- ls(pattern = "dataset_",envir = environment())
sim_function = sim_custom2; lik_function = lik_custom
Nsims <- 1E6
which_trees <- 1:length(data.sets); #which_trees <- c(14); i <- 1 #use it in case you already have some partial results
for (i in which_trees)
{
  print(i)
  d.s <- get(data.sets[[i]])
  #the aim is to get "lik_result" equal to "sim_result" for an high enough number of simulations
  test_result1 <- test_likelihood_formula(dataset = d.s, Nsims = Nsims, lik_function = lik_function, sim_function = sim_function); print(test_result1$results.table)
  if (Nsims >= 1E5)
  {
    results_file <- paste0(getwd(),"//results//table3.xls")
    sheet_name <- sheet_name0 <- data.sets[[i]]
    wb <- xlsx::loadWorkbook(file = results_file)
    prova <- xlsx::getSheets(wb = wb)
    nomi <- names(prova)
    jj <- 1; while (sheet_name %in% nomi) {jj <- jj + 1; sheet_name <- paste0(sheet_name0," - ", toString(jj))}

    xlsx::write.xlsx(x = test_result1$results.table, file = results_file, sheetName = sheet_name, append = TRUE)
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














