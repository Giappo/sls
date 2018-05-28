rm(list = ls())
library(sls)
source(paste0(getwd(),"//R//datasets.R"))
sim_function = sim_custom; lik_function = lik_custom
Nsims <- 1000000 # Nsims <- 10000000

d.s <- dataset_pure_branching3
#the aim is to get "lik_result" equal to "sim_result" for an high enough number of simulations
test_result1 <- test_likelihood_formula2(dataset = d.s, Nsims = Nsims, lik_function = lik_function, sim_function = sim_function); print(test_result1$results.table)
if (Nsims >= 100000 && all.equal(sim_function, sim_custom) && all.equal(lik_function, lik_custom))
{
  results_file <- paste0(getwd(),"//results//table.xls")
  xlsx::write.xlsx(x = test_result1$results.table, file = results_file, sheetName = test_result1$sheet_name, append = TRUE)
  # xlsx::addPicture(file = results_file, sheet = test_result$sheet_name, startRow = 13, startColumn = 1)
}

Nsims <- 1000000 # Nsims <- 10000000
sim_function = sim_R_example; lik_function = lik_custom
d.s <- dataset_Rampal
#the aim is to get "lik_result" equal to "sim_result" for an high enough number of simulations
test_result2 <- test_likelihood_formula2(dataset = d.s, Nsims = Nsims, lik_function = lik_function, sim_function = sim_function); print(test_result2$results.table)
if (Nsims >= 100000 && all.equal(sim_function, sim_custom) && all.equal(lik_function, lik_custom))
{
  results_file <- paste0(getwd(),"//results//table.xls")
  xlsx::write.xlsx(x = test_result2$results.table, file = results_file, sheetName = test_result2$sheet_name, append = TRUE)
  # xlsx::addPicture(file = results_file, sheet = test_result$sheet_name, startRow = 13, startColumn = 1)
}

# file <- system.file("tests", "log_plot.jpeg", package = "xlsx")
# file <- "accuracy_vs_samplesize.png"
# # wb <- xlsx::createWorkbook()
# wb <- results_file
# xlsx::getSheets(wb)
# sheet <- createSheet(wb, test_result$sheet_name)
# addPicture(file, sheet)
