# rm(list = ls())
# functions
percentiles_function <- function(results, sim_pars, printit = 1,
                                 quantiles_choice = c(.25, .50, .75)){
  quantiles_names <- format(round(quantiles_choice, 2), nsmall = 2)
  Npars <- length(sim_pars)
  parnames <- colnames(results)[1:Npars]
  percentiles <- vector("list", 3)
  for (idpar in 1:Npars){percentiles[[idpar]] <- stats::quantile(results[,idpar], quantiles_choice)}
  percentiles <- t(matrix(unlist(percentiles), nrow = 3, ncol = Npars));
  colnames(percentiles) <- quantiles_names; rownames(percentiles) <- parnames
  if (printit==1)
  {
    print(percentiles); print(sim_pars)
  }
  out <- percentiles
}

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

model_comparison <- function(sls.Results, DDD.Results, dataset_pars1, path) { #title,

  multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL, title = "",
                        fontsize = 12, fontfamily = "Helvetica") {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (nchar(title)>0){
      layout <- rbind(rep(0, ncol(layout)), layout)
    }

    if (numPlots==1) {
      print(plots[[1]])

    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout),
                                                 ncol(layout),
                                                 heights = if (nchar(title)>0) {unit(c(0.5, rep(5,nrow(layout)-1)), "null")}
                                                 else {unit(c(rep(5, nrow(layout))), "null")})))

      # Make each plot, in the correct location
      if (nchar(title)>0) {
        grid.text(title,
                  vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol(layout)),
                  gp = gpar(fontsize = fontsize, fontfamily = fontfamily))
      }

      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }

df <- data.frame(sls.Results, DDD.Results)
names(df) <- c(paste0("sls.", colnames(sls.Results)), paste0("DDD.", colnames(DDD.Results)))

pars_string <- gsub(", ", "_", toString(dataset_pars1))

pdfname <- paste0("models_comparison_[lambda-mu]_", pars_string)
grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
p5 <- ggplot2::ggplot(df,
                      ggplot2::aes(x = (sls.lambda_M - sls.mu_M),
                                   y = (DDD.lambda_M - DDD.mu_M)
                                   )
      ) +
      ggplot2::geom_point(color = "firebrick4"
      ) +
      ggplot2::labs(x = expression(lambda[M]^sls - mu[M]^sls),
                    y = expression(lambda[M]^DDD - mu[M]^DDD)
      ) +
      ggplot2::theme_bw(
      ); print(p5)
grDevices::dev.off()

pdfname <- paste0("models_comparison_[mulambda_ratio]_", pars_string)
grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
p6 <- ggplot2::ggplot(df,
                      ggplot2::aes(x = (sls.mu_M / sls.lambda_M),
                                   y = (DDD.mu_M / DDD.lambda_M)
                      )
      ) +
      ggplot2::geom_point(color = "turquoise4"
      ) +
      ggplot2::labs(x =  expression(mu[M]^sls/lambda[M]^sls) ,
                    y =  expression(mu[M]^DDD/lambda[M]^DDD)
      ) +
      ggplot2::theme_bw(
      ); print(p6)
grDevices::dev.off()

pdfname <- paste0("models_comparison_[lambda]_", pars_string)
grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
p7 <- ggplot2::ggplot(df,
                      ggplot2::aes(x = (sls.lambda_M),
                                   y = (DDD.lambda_M)
                      )
      ) +
      ggplot2::geom_point(color = "firebrick4"
      ) +
      ggplot2::labs(x = expression(lambda[M]^sls),
                    y = expression(lambda[M]^DDD)
      ) +
      ggplot2::theme_bw(
      ); print(p7)
grDevices::dev.off()

pdfname <- paste0("models_comparison_[mu]_", pars_string)
grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
p8 <- ggplot2::ggplot(df,
                      ggplot2::aes(x = (sls.mu_M),
                                   y = (DDD.mu_M)
                      )
      ) +
      ggplot2::geom_point(color = "turquoise4"
      ) +
      ggplot2::labs(x =  expression(mu[M]^sls) ,
                    y =  expression(mu[M]^DDD)
      ) +
      ggplot2::theme_bw(
      ); print(p8)
grDevices::dev.off()

# pp <- multiplot(p1, p2, p3, p4, cols = 2, title = title)
# pp <- multiplot(p5, p6, cols = 2, title = title)

return()
}

showplot.box <- function(results) {}

showplot.correlation <- function(results) {}

showplot.cloud <- function(results) {}

sls.analyze_data <- function(max_sims = 1000,
                             quantiles_choice = c(.25, .50, .75),
                             results_folder = "results"){

  home_dir <- substring(getwd(), 1, 21)
  proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir, "/Progress/")))
  folder_name <- paste0(home_dir, "/Progress/", list.files(paste0(home_dir, "/Progress"))[proj.coords])
  # project_folder <- paste0(folder_name, "/sls")
  results_folder <- paste0(folder_name, "/results")
  datasets <- list.files(results_folder, pattern = "^[0]")

  dataset_pars <- vector("list", Nd <- length(datasets))
  for (d in 1:Nd)
  {
    dataset_pars[[d]] <- as.numeric(unlist( strsplit(x = datasets[d], split = "-")))
  }

  parnames = c("lambda_M", "mu_M", "K_M", "lambda_S", "mu_S", "K_S", "t_d")
  idparsopt <- c(1,2,4,5); Npars <- length(idparsopt)
  sls_quantiles <- DDD_quantiles <- quantiles <- vector("list", Nd); Nsims <- rep(0, Nd)
  for (d in 1:Nd)
  {
    print(d)
    sim_pars <- dataset_pars[[d]][1:Npars];
    cond <- dataset_pars[[d]][Npars + 1]
    local_path <- paste0(results_folder,"//", datasets[d]);

    DDD_files <- list.files(pattern = paste0('DDD_MLE'), path = local_path, full.names = TRUE); length(DDD_files)
    sls_files <- list.files(pattern = paste0('sls_MLE'), path = local_path, full.names = TRUE); length(sls_files)

    if (length(sls_files) != 0 && length(DDD_files) != 0)
    {
      for (s in 1:length(sls_files))
      {
        fileData <- read.table(file = sls_files[s], header = FALSE, sep = ",")
        ifelse(exists("targetTable"), targetTable <- rbind(targetTable, fileData), targetTable <- fileData)
      }
      sls_results0 <- targetTable; rm(targetTable); names(sls_results0) <- (c(parnames, "LL", "df", "conv", "tree_id"))

      for (s in 1:length(DDD_files))
      {
        fileData <- read.table(file=DDD_files[s], header = FALSE, sep = ",")
        ifelse(exists("targetTable"),targetTable <- rbind(targetTable, fileData), targetTable <- fileData)
      }
      DDD_results0 <- targetTable; rm(targetTable); names(DDD_results0) <- (c(parnames, "LL", "df", "conv", "tree_id"))

      sls_right_ids  <- sls_results0[rowSums(sls_results0[, idparsopt] == rep(-1, Npars)) != Npars, "tree_id"]; length(sls_right_ids)
      DDD_right_ids  <- DDD_results0[rowSums(DDD_results0[, idparsopt] == rep(-1, Npars)) != Npars, "tree_id"]; length(DDD_right_ids)
      right_ids <- intersect(sls_right_ids, DDD_right_ids); length(right_ids)

      sls_results1 <- sls_results0[sls_results0[, "tree_id"] %in% right_ids,]; dim(sls_results1)
      sls_results  <- sls_results1[,idparsopt]; Nsims[d] <- nrow(sls_results)
      if (max(sls_results[, 1:Npars] == -1)){print("You are considering results that are = -1. Be careful!")}

      sls_quantiles[[d]] <- percentiles_function(results = sls_results, sim_pars = sim_pars,
                                                       printit = 0, quantiles_choice = quantiles_choice)
      titolo  <- paste0("sls - Correlation analysis with ", Nsims[d], "/", max_sims, " trees.")
      pdfname <- paste0("sls_Correlation ", datasets[d])
      correlation_analysis(results = sls_results, sim_pars = sim_pars, titolo = titolo,
                           pdfname = pdfname, path = results_folder, openit = 0)

      DDD_results1 <- DDD_results0[DDD_results0[, "tree_id"] %in% right_ids,]; dim(DDD_results1)
      DDD_results  <- DDD_results1[,idparsopt]; Nsims[d] <- nrow(DDD_results)
      if (max(DDD_results[, 1:Npars] == -1) ){print("You are considering results that are = -1. Be careful!")}

      DDD_quantiles[[d]] <- percentiles_function(results = DDD_results, sim_pars = sim_pars,
                                                       printit = 0, quantiles_choice = quantiles_choice)
      titolo  <- paste0("DDD - Correlation analysis with ", Nsims[d], "/", max_sims, " trees.")
      pdfname <- paste0("DDD_Correlation ", datasets[d])
      correlation_analysis(results = DDD_results, sim_pars = sim_pars, titolo = titolo,
                           pdfname = pdfname, path = results_folder, openit = 0)

      model_comparison(sls.Results = sls_results1,
                       DDD.Results = DDD_results1,
                       dataset_pars1 = dataset_pars[[d]],
                       path = results_folder)
    }
  }

  parnames2 <- parnames[idparsopt]
  sls_result.table <- matrix(NA, nrow = Nd, ncol = 1 + prod(dim(sls_quantiles[[1]])) + Npars)
  DDD_result.table <- matrix(NA, nrow = Nd, ncol = 1 + prod(dim(DDD_quantiles[[1]])) + Npars)
  quantiles_names  <- format(round(quantiles_choice, 2), nsmall = 2)

  par_quantiles_names <- vector("list", Npars)
  for (ll in seq_along(par_quantiles_names))
  {
    # par_quantiles_names[[ll]] <- paste0(parnames2[ll], quantiles_names)
    par_quantiles_names[[ll]] <- paste0(parnames2[ll], quantiles_names)
  }; par_quantiles_names

  for (d in 1:(Nd))
  {
    # sls_result.table[d, 1:Npars] <- dataset_pars[[d]][1:Npars]
    sls_result.table[d, 1:(Npars + 1)] <- dataset_pars[[d]][1:(Npars + 1)]
    if (Nsims[d] != 0)
    {
      # sls_result.table[d, (Npars + 1):(Npars + prod(dim(sls_quantiles[[1]])))] <-
      sls_result.table[d, (1 + Npars + 1):(1 + Npars + prod(dim(sls_quantiles[[1]])))] <-
        t(matrix(t(sls_quantiles[[d]]), nrow = prod(dim(sls_quantiles[[1]])), byrow = F));
    }
  }
  for (d in 1:(Nd))
  {
    # DDD_result.table[d, 1:Npars] <- dataset_pars[[d]][1:Npars]
    DDD_result.table[d, 1:(Npars + 1)] <- dataset_pars[[d]][1:(Npars + 1)]
    if (Nsims[d] != 0)
    {
      # DDD_result.table[d, (Npars + 1):(Npars + prod(dim(DDD_quantiles[[1]])))] <-
      DDD_result.table[d, (1 + Npars + 1):(1 + Npars + prod(dim(DDD_quantiles[[1]])))] <-
        t(matrix(t(DDD_quantiles[[d]]), nrow = prod(dim(DDD_quantiles[[1]])), byrow = F));
    }
  }

  sim_par_names <- paste0("sim.", parnames2)
  # colnames(sls_result.table) <- colnames(DDD_result.table) <- c(sim_par_names, unlist(par_quantiles_names))
  colnames(sls_result.table) <- colnames(DDD_result.table) <- c(sim_par_names, "cond", unlist(par_quantiles_names))
  sls_result.table2 <- format(round(sls_result.table, 2), nsmall = 2)
  sls_result.table2 <- cbind(sls_result.table2, Nsims); sls_result.table2
  DDD_result.table2 <- format(round(DDD_result.table, 2), nsmall = 2)
  DDD_result.table2 <- cbind(DDD_result.table2, Nsims); DDD_result.table2
  # install.packages("xlsx");
  library("xlsx")
  xlsname0 = "results_table"; xlsname = xlsname0; xlscount = 2;
  while (file.exists(paste0(folder_name,"//",xlsname,".xlsx"))){xlsname = paste0(xlsname0,xlscount); xlscount = xlscount + 1}
  write.xlsx(x = sls_result.table2, file = paste0(folder_name,"//", xlsname,".xlsx"), sheetName = paste0("sls - age=10; #sims=", max_sims), row.names = FALSE)
  write.xlsx(x = DDD_result.table2, file = paste0(folder_name,"//", xlsname,".xlsx"), sheetName = paste0("DDD - age=10; #sims=", max_sims), row.names = FALSE, append = TRUE)
  return(list(sls_result.table = sls_result.table, DDD_result.table = DDD_result.table))
}

sls.analyze_data2 <- function(resultstable, error_bars = 0){

  library(ggplot2); library(RColorBrewer)

  pippo_draw <- function(sub_res_table, error_bars = error_bars, model_name){

    sub_res_table <- data.frame(sub_res_table)

    if (1) { #2nd attempt
      library(ggplot2)

      data2 <- sub_res_table
      Npars <- 4
      true_values <- unique(data2[,1:Npars])
      data3 <- data2[1:nrow(true_values), 1:ncol(data2)]
      for (i in 1:Npars)
      {
        data3[,(3:5) + 3*i] <- true_values[,i]
      }
      data3[,1:Npars] <- true_values
      data3$cond <- -1
      data4 <- rbind(data2, data3)

      colors <- factor(data4$sim.lambda_M^0.98 * data4$sim.mu_M^1.02);
      colors_names <- apply(cbind(data4$sim.lambda_M, data4$sim.mu_M), MARGIN = 1, FUN = "toString")
      # xx <- paste0(expression(lambda[M]), "=", data4$sim.lambda_M); yy <- paste0(expression(mu[M]), "=", data4$sim.mu_M)
      # # xx <- paste0(expression("λ[M]"), " = ", data4$sim.lambda_M); # yy <- paste0(expression("μ[M]"), " = ", data4$sim.mu_M)
      # colors_names <- apply(cbind(xx, yy), MARGIN = 1, FUN = function(x) paste0(x[1], ", ", x[2]))

      palette  <- rainbow(n = nrow(sub_res_table), s = 0.5)
      palette2 <- rainbow(n = length(unique(colors)), s = 0.5)

      text_size <- 12
      sp4 <- ggplot(data4,
                    aes(x = lambda_M0.50,
                        y = mu_M0.50,
                        shape = factor(data4$cond),
                        col = colors,
                        size = 0.05
                    )
             ) +
             ggtitle(paste0("Results for model ", model_name)
             ) +
             ggplot2::labs(x = expression(lambda[M]), #"lambda_M"
                           y = expression(mu[M])#"mu_M"
             ) +
             geom_point(
             ) +
             scale_color_manual(name   = bquote("Parameter setting \n(" ~lambda[M]~","~mu[M]~")"),
                                values = palette2,
                                labels = unique((colors_names))
             ) +
             scale_shape_manual(name   = "Conditioning",
                                values = c(19, 21, 24, 22),
                                labels = c("True values", "0", "1", "2")
             ) +
             guides(size = FALSE #remove "size" from the legend
             ) +
             theme(legend.text  = element_text(size = (text_size)),
                   legend.title = element_text(size = (text_size + 4), face = "bold"),
                   axis.title   = element_text(size = (text_size + 2), face = "bold"),
                   axis.text    = element_text(size = text_size),
                   plot.title   = element_text(size = 18, face = "bold")
             ) +
             theme_dark(
             ) +
             guides(shape = guide_legend(override.aes = list(size = 3))); sp4# scatterplot
    }

    return(sp4)
  }

  Lr <- length(result.table)
  model_names <- c("sls", "DDD"); ii <- 1

  for (ii in 1:Lr)
  {
    res <- resultstable[[ii]][apply(!is.na(resultstable[[ii]]), 1, prod) != 0,]
    pippo <- as.data.frame(res)
    colnames(pippo) <- colnames(res)
    pippo2 <- cbind(pippo, row(pippo)[,1])

    ###
    folder_name <- "F://Dropbox//University//Progress//RQ4-single_lineage_rate_shifts//results"
    setwd(folder_name)

    # grDevices::png(filename = paste0(model_names[[ii]], "_results.png"))
    grDevices::pdf(file = paste0(model_names[[ii]], "_results.pdf"))
    results_plot <- pippo_draw(sub_res_table = pippo2, model_name = model_names[[ii]])
    print(results_plot)
    grDevices::dev.off()
  }
}

# actual analysis
resultstable <- sls.analyze_data(max_sims = 2000)
sls.analyze_data2(resultstable = resultstable)
