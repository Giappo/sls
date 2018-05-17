rm(list = ls())
library(sls)
source(paste0(getwd(),"//R//datasets.R"))

Nsims <- 10000000 # Nsims <- 10000000
d.s <- dataset_1
# d.s <- dataset_2
# d.s <- dataset_Rampal
# d.s <- dataset_Bart

#the aim is to get "lik_result" equal to "sim_result" for an high enough number of simulations
test_result <- test_likelihood_formula2(dataset = d.s, Nsims = Nsims); print(test_result$results.table)
if (Nsims >= 100000)
{
  results_file <- paste0(getwd(),"//results//table.xls")
  xlsx:::write.xlsx(x = test_result$results.table, file = results_file, sheetName = test_result$sheet_name, append = TRUE)
  xlsx:::addPicture(file = results_file, sheet = test_result$sheet_name, startRow = 13, startColumn = 1)
}

file <- system.file("tests", "log_plot.jpeg", package = "xlsx")
file <- "accuracy_vs_samplesize.png"

# wb <- xlsx:::createWorkbook()
wb <- results_file
xlsx:::getSheets(wb)
sheet <- createSheet(wb, test_result$sheet_name)

addPicture(file, sheet)
