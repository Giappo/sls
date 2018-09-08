correlation_analysis <- function(results, path, titolo = NULL, pdfname, sim_pars = sim_pars,
         percentage_hidden_outliers = 0, openit = 0){

  Npars <- length(sim_pars);
  par_names <- colnames(results)

  truevalues_color = "red"; points_color = "azure4"; medians_color = "blue3"#"chartreuse3";
  medians_color_name = "blue"; truevalues_color_name = "red";

  medians <- rep(0, Npars); for (idpar in 1:Npars){medians[idpar] <- median(results[, idpar])}
  medians_string <- paste0( "MLE Medians (",medians_color_name,") = (",paste(signif(medians, 2),sep = "''",collapse = ", "),")")
  truevalues_string <- paste0( "True Values (",truevalues_color_name,") = (",paste(signif(sim_pars,2),sep = "''",collapse = ", "),")")
  axislimits <- rep(NA, Npars)
  for (i in 1:Npars)
  {
    axislimits[i] <- quantile(results[,i], probs = 1 - percentage_hidden_outliers)
  }

  #pdf
  pdf(file = paste0(path, "/", pdfname, ".pdf"));
  par(mfrow = c(Npars, Npars)); par(oma = c(0,0,2,0));
  for (i in 1:Npars){for (j in 1:Npars){
    good.lines = results[,i] < axislimits[i] & results[,j] < axislimits[j]
    ifelse(any(good.lines) > 0,good.results <- results[good.lines,], good.results <- results)

    if (i == j){hist((good.results[,i]), main = NULL, xlab = paste(par_names[i]), breaks = 15); #,breaks = 15
      abline(v = sim_pars[i], col = truevalues_color)
      abline(v = medians[i], col = medians_color)
    }
    else{plot(good.results[,i] ~ good.results[,j],xlab=par_names[j],ylab=par_names[i],cex=0.3,col=points_color);
      points(x = sim_pars[j], y = sim_pars[i], col = truevalues_color, pch = 10, cex = 1.5)
      points(x = medians[j] , y = medians[i] , col = medians_color, pch = 10, cex = 1.5)
    }
  }}
  title(main = (titolo.pdf <- (paste0("\n\n", titolo, "\n", medians_string, "\n", truevalues_string))), outer = TRUE)
  dev.off()

  if (openit == 1)
  {
    file.show(normalizePath(paste0(path, "/", pdfname, ".pdf")))
  }
}


sls.analyze_data <- function(max_sims = 1000,
                             quantiles_choice = c(.25, .50, .75),
                             results_folder = "results"){

  home_dir <- substring(getwd(), 1, 21)
  proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir,"/Progress/")))
  folder_name <- paste0(home_dir,"/Progress/", list.files(paste0(home_dir,"/Progress"))[proj.coords])
  project_folder <- paste0(folder_name, "/sls")
  results_folder <- paste0(project_folder, "/results")
  datasets <- list.files(results_folder, pattern = "^[0]")

  dataset_pars <- vector("list",Nd <- length(datasets))
  for (d in 1:Nd)
  {
    dataset_pars[[d]] <- as.numeric(  unlist( strsplit(x = datasets[d],split = "-") )  )
  }

  parnames = c("lambda_M", "mu_M", "K_M", "lambda_S", "mu_S", "K_S", "t_d")
  idparsopt <- c(1,2,4,5); Npars <- length(idparsopt)
  quantiles <- vector("list", Nd); Nsims <- rep(0, Nd)
  for (d in 1:Nd){
    print(d)
    sim_pars = dataset_pars[[d]][1:Npars];
    cond = dataset_pars[[d]][Npars + 1]
    local_path <- paste0(results_folder,"//", datasets[d]);

    DDD_files <- list.files(pattern = paste0('DDD_MLE'), path = local_path, full.names = TRUE); length(DDD_files)
    sls_files <- list.files(pattern = paste0('sls_MLE'), path = local_path, full.names = TRUE); length(sls_files)

    res_files <- sls_files
    if (length(res_files) != 0)
    {
      for (s in 1:length(res_files))
      {
        fileData <- read.table(file = res_files[s], header = FALSE, sep = ",")
        ifelse(exists("targetTable"),targetTable <- rbind(targetTable, fileData), targetTable <- fileData)
      }

      results0 <- targetTable; rm(targetTable); names(results0) <- (c(parnames, "LL", "df", "conv", "tree_id"))
      results1  <- results0[rowSums(results0[, idparsopt] == rep(-1, Npars)) != Npars,];
      results   <- results1[,idparsopt]

      Nsims[d] <- nrow(results)
      if ( max(results[, 1:Npars] == -1) ){print("You are considering results that are = -1. Be careful!")}

      quantiles[[d]] <- MBD:::percentiles_function(results = results, sim_pars = sim_pars,
                                                   printit = 0, quantiles_choice = quantiles_choice)
      titolo <- paste("sls - Correlation analysis with ", Nsims[d],"/",max_sims," trees.",sep = '')
      pdfname <- paste("sls_Correlation ", datasets[d], sep = '')
      correlation_analysis(results = results, sim_pars = sim_pars, titolo = titolo,
                                 pdfname = pdfname, path = results_folder, openit = 0)
    }

    res_files <- DDD_files
    if(length(res_files) != 0)
    {
      for (s in 1:length(res_files))
      {
        fileData <- read.table(file=res_files[s], header = FALSE, sep = ",")
        ifelse(exists("targetTable"),targetTable <- rbind(targetTable, fileData), targetTable <- fileData)
      }

      results0 <- targetTable; rm(targetTable); names(results0) <- (c(parnames, "LL", "df", "conv", "tree_id"))
      results1  <- results0[rowSums(results0[, idparsopt] == rep(-1, Npars)) != Npars,];
      results   <- results1[,idparsopt]

      Nsims[d] <- nrow(results)
      if ( max(results[, 1:Npars] == -1) ){print("You are considering results that are = -1. Be careful!")}

      quantiles[[d]] <- MBD:::percentiles_function(results = results, sim_pars = sim_pars,
                                                   printit = 0, quantiles_choice = quantiles_choice)
      titolo <- paste("DDD - Correlation analysis with ", Nsims[d],"/",max_sims," trees.",sep = '')
      pdfname <- paste("DDD_Correlation ", datasets[d], sep = '')
      correlation_analysis(results = results, sim_pars = sim_pars, titolo = titolo,
                                 pdfname = pdfname, path = results_folder, openit = 0)
    }
  }

  parnames2 <- parnames[idparsopt]
  result.table <- matrix(NA,nrow = Nd,ncol = prod(dim(quantiles[[1]])) + Npars )
  quantiles_names <- format(round(quantiles_choice,2),nsmall = 2)

  par_quantiles_names <- vector("list", Npars)
  for (ll in seq_along(par_quantiles_names))
  {
    par_quantiles_names[[ll]] <- paste0(parnames2[ll], quantiles_names)
  }

  for (d in 1:(Nd))
  {
    result.table[d, 1:Npars] <- dataset_pars[[d]][1:Npars]
    if (Nsims[d] != 0)
    {
      result.table[d, (Npars + 1):(Npars + prod(dim(quantiles[[1]])))] = t(matrix(t(quantiles[[d]]), nrow = prod(dim(quantiles[[1]])), byrow = F));
    }
  }
  sim_par_names <- paste0("sim.", parnames2)
  colnames(result.table) = c(sim_par_names, unlist(par_quantiles_names))
  result.table2 <- format(round(result.table,2),nsmall = 2)
  result.table2 <- cbind(result.table2,Nsims); result.table2
  # install.packages("xlsx");
  library("xlsx")
  xlsname0 = "results_table"; xlsname = xlsname0; xlscount = 2;
  while (file.exists(paste0(folder_name,"//",xlsname,".xlsx"))){xlsname = paste0(xlsname0,xlscount); xlscount = xlscount + 1}
  write.xlsx(x = result.table2, file = paste(folder_name,"//", xlsname,".xlsx",sep = ''), sheetName = paste0("age=10; #sims=", max_sims), row.names = FALSE)
  return(result.table)
}

sls.analyze_data(max_sims = 1500)
