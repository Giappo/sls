# rm(list = ls())
# functions
collect_data <- function() {

  get.dropbox.folder <- function() {

    list.of.packages <- c("RJSONIO")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
    require("RJSONIO")

    if (Sys.info()['sysname'] == 'Darwin')
    {
      info <- RJSONIO::fromJSON(
        file.path(path.expand("~"),'.dropbox','info.json'))
    }
    if (Sys.info()['sysname'] == 'Windows')
    {
      info <- RJSONIO::fromJSON(
        if (file.exists(file.path(Sys.getenv('APPDATA'), 'Dropbox','info.json'))) {
          file.path(Sys.getenv('APPDATA'), 'Dropbox', 'info.json')
        } else {
          file.path(Sys.getenv('LOCALAPPDATA'),'Dropbox','info.json')
        }
      )
    }
    dropbox_base <- info$personal$path
  }
  db_dir <- get.dropbox.folder()
  home_dir <- db_dir
  home_dir <- paste0(db_dir, "\\university\\Progress\\"); list.files(home_dir)
  proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir)))
  folder_name <- paste0(home_dir, list.files(paste0(home_dir))[proj.coords])
  results_mother_folder <- paste0(folder_name, "/results")
  results_folder <- choose.dir(default = results_mother_folder)
  datasets <- list.files(results_folder, pattern = "^[0]")

  dataset_pars <- vector("list", Nd <- length(datasets))
  for (d in 1:Nd)
  {
    dataset_pars[[d]] <- as.numeric(unlist( strsplit(x = datasets[d], split = "-")))
  }

  parnames  <- c("lambda_M", "mu_M", "K_M", "lambda_S", "mu_S", "K_S", "t_d")
  useful_quantities <- c("lambda_M", "mu_M", "lambda_S", "mu_S", "LL", "tree_id")
  idparsopt <- c(1,2,4,5); Npars <- length(idparsopt)
  parnames2 <- parnames[idparsopt]
  sim.parnames <- paste0("sim.", parnames2)
  MLE.parnames <- paste0("MLE.", parnames2)
  all_results  <- NULL; d <- 1
  for (d in 1:Nd)
  {# d loop
    print(d)

    sim_pars   <- dataset_pars[[d]][1:Npars];
    cond       <- dataset_pars[[d]][Npars + 1]
    local_path <- paste0(results_folder, "//", datasets[d]);
    names0  <- list.files(pattern = "[.]txt", path = local_path, full.names = TRUE)
    names1 <- gsub(".*/","", names0)
    model_names <- unique(gsub("_MLE.*","", names1))
    model_ids  <- model_results2 <- model_results <- vector("list", length(model_names)); m <- 1
    d_results  <- NULL; m <- 1
    for (m in seq_along(model_names))
    {# m loop
      filenames <- paste0(model_names[m], "_MLE")
      files <- list.files(pattern = paste0(filenames), path = local_path, full.names = TRUE); length(files)

      if (length(files) != 0)
      {
        for (s in 1:length(files))
        {# s loop
          fileData <- read.table(file = files[s], header = FALSE, sep = ",")
          ifelse(exists("targetTable"), targetTable <- rbind(targetTable, fileData), targetTable <- fileData)
        }# s loop
        results0  <- targetTable; rm(targetTable);
        names(results0) <- c(parnames, "LL", "df", "conv", "tips_M", "tips_S", "tree_id")[1:ncol(results0)]
        names(results0)[ncol(results0)] <- "tree_id"

        results1  <- results0[, names(results0) %in% useful_quantities]; names(results1)
        results2  <- results1; names(results2)[1:Npars] <- MLE.parnames
        right_ids <- results2[results2[, 1] != -1, "tree_id"]; length(right_ids)
        results3  <- results2[results2[, "tree_id"] %in% right_ids,]; dim(results3)
        results4  <- results3[order(results3$tree_id),]
        results5  <- cbind(results4,
                           rep(model_names[m], nrow(results4)),
                           matrix(dataset_pars[[d]],
                                  ncol = length(dataset_pars[[d]]),
                                  nrow = nrow(results4),
                                  byrow = TRUE)
                           )
        names(results5) <- c(names(results4), "model", sim.parnames, "cond"); head(results5)
        testit::assert(results5$cond == floor(results5$cond) | results5$cond == 0)
      }
      model_ids[[m]] <- right_ids
      d_results <- rbind(d_results, results5)
      rm(results0, results1, results2, results3, results4, results5,
         right_ids, files, filenames)
    }# m loop
    d_ids <- Reduce(intersect, model_ids); d_ids <- sort(d_ids)
    d_results <- d_results[(d_results[,"tree_id"] %in% d_ids),]
    all_results <- rbind(all_results, d_results)
    rm(d_results)

  }# d loop

  for (m in seq_along(model_names)[-1])
  {
    testit::assert(
    dim(all_results[all_results$model == model_names[m - 1],]) == dim(all_results[all_results$model == model_names[m],])
    )
  }

  all_results$directory <- results_folder

  return(all_results)
}

showplot.box <- function(data) {

  library(ggplot2); library(reshape2); library(scales)

  # home_dir <- substring(getwd(), 1, 21)
  # proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir, "/Progress/")))
  # folder_name <- paste0(home_dir, "/Progress/", list.files(paste0(home_dir, "/Progress"))[proj.coords])
  # results_folder <- paste0(folder_name, "/results")
  results_folder <- data$directory[1]
  path <- paste0(results_folder, "/boxplots")
  if (!file.exists(path)) {dir.create(file.path(path), showWarnings = FALSE)}

  #   geom_boxplot(outlier.colour = "black", outlier.shape = 16,
  #                outlier.size = 2, notch = FALSE)

  data1 <- data
  data1$sim.lamu <- interaction(data1$sim.lambda_M, data1$sim.mu_M)
  data1.m <- melt(data1, id.vars = 'sim.lamu',
                  measure.vars = "model")
  data1.m$delta_lambda_M  <- data1$MLE.lambda_M - data1$sim.lambda_M
  data1.m$delta_lambda_M2 <- (data1$MLE.lambda_M - data1$sim.lambda_M)/data1$sim.lambda_M
  data1.m$delta_mu_M      <- data1$MLE.mu_M - data1$sim.mu_M
  data1.m$delta_mu_M2     <- (data1$MLE.mu_M - data1$sim.mu_M)/data1$sim.mu_M
  data1.m$delta_mulambda  <- (data1$MLE.mu_M/data1$MLE.lambda_M) - (data1$sim.mu_M/data1$sim.lambda_M)
  data1.m$model <- data1$model
  data1.m$cond  <- data1$cond
  # levels(data1.m$cond) <- c("No Cond", "Cond 1", "Cond 2")
  xlabels <- unique(paste0(data1$sim.lambda_M, "\n", data1$sim.mu_M))

  pl <- ggplot(data1.m
  ) +
  geom_boxplot(aes(x = sim.lamu, y = delta_lambda_M, color = model)
    # geom_boxplot(aes(x = sim.lamu, y = delta_lambda_M2, color = model)
  ) +
  facet_grid(. ~ cond
  ) +
  xlab(bquote("Parameter setting" ~ "(" ~ lambda[M] ~ "," ~ mu[M] ~ ")")
  ) +
  ylab(bquote(Delta ~ lambda[M])
  ) +
  theme(axis.line = element_line(colour = "darkblue",
                                 size = 1, linetype = "solid")
  ) +
  scale_x_discrete(labels = xlabels); pl

  pm <- ggplot(data1.m
  ) +
  geom_boxplot(aes(x = sim.lamu, y = delta_mu_M, color = model)
    # geom_boxplot(aes(x = sim.lamu, y = delta_mu_M2, color = model)
  ) +
  facet_grid(. ~ cond
  ) +
  xlab(bquote("Parameter setting" ~ "(" ~ lambda[M] ~ "," ~ mu[M] ~ ")")
  ) +
  ylab(bquote(Delta ~ mu[M])
  ) +
  theme(axis.line = element_line(colour = "darkblue",
                                 size = 1, linetype = "solid")
  ) +
  scale_x_discrete(labels = xlabels); pm

  plm1 <- ggplot(data1.m
  ) +
  geom_boxplot(aes(x = sim.lamu, y = (delta_lambda_M - delta_mu_M), color = model)
               # geom_boxplot(aes(x = sim.lamu, y = delta_mu_M2, color = model)
  ) +
  facet_grid(. ~ cond
  ) +
  xlab(bquote("Parameter setting" ~ "(" ~ lambda[M] ~ "," ~ mu[M] ~ ")")
  ) +
  ylab(bquote(Delta ~ (lambda[M] - mu[M]))
  ) +
  theme(axis.line = element_line(colour = "darkblue",
                                 size = 1, linetype = "solid")
  ) +
  scale_x_discrete(labels = xlabels); plm1

  plm2 <- ggplot(data1.m
  ) +
  geom_boxplot(aes(x = sim.lamu, y = delta_mulambda, color = model)
               # geom_boxplot(aes(x = sim.lamu, y = delta_mu_M2, color = model)
  ) +
  facet_grid(. ~ cond
  ) +
  xlab(bquote("Parameter setting" ~ "(" ~ lambda[M] ~ "," ~ mu[M] ~ ")")
  ) +
  ylab(bquote(Delta ~ (mu[M]/lambda[M]))
  ) +
  theme(axis.line = element_line(colour = "darkblue",
                                 size = 1, linetype = "solid")
  ) +
  scale_x_discrete(labels = xlabels); plm2

  pdfname <- paste0("boxplot_lambdaM")
  grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
  print(pl)
  grDevices::dev.off()

  pdfname <- paste0("boxplot_muM")
  grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
  print(pm)
  grDevices::dev.off()

  pdfname <- paste0("boxplot_lambdaM-muM")
  grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
  print(plm1)
  grDevices::dev.off()

  pdfname <- paste0("boxplot_ratiomuMlambdaM")
  grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
  print(plm2)
  grDevices::dev.off()

  return(list(pl, pm, plm1, plm2))
}

showplot.correlation <- function(data) {

  library(ggplot2); library(reshape2); library(scales)

  # home_dir <- substring(getwd(), 1, 21)
  # proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir, "/Progress/")))
  # folder_name <- paste0(home_dir, "/Progress/", list.files(paste0(home_dir, "/Progress"))[proj.coords])
  # results_folder <- paste0(folder_name, "/results")
  results_folder <- data$directory[1]
  path <- paste0(results_folder, "/correlationplots")
  if (!file.exists(path)) {dir.create(file.path(path), showWarnings = FALSE)}

  data2 <- data
  data2$parsetting <- interaction(data2$sim.lambda_M, data2$sim.mu_M,
                                  data2$sim.lambda_S, data2$sim.mu_S,
                                  data2$cond)

  ds <- data2[data2$model == "sls",]; names(ds) <- paste0("sls.", names(data2))
  dD <- data2[data2$model == "DDD",]; names(dD) <- paste0("DDD.", names(data2))
  df <- data.frame(ds, dD); names(df)
  df$parsetting <- interaction(df$sls.sim.lambda_M, df$sls.sim.mu_M,
                               df$sls.sim.lambda_S, df$sls.sim.mu_S,
                               df$sls.cond)

  for (var in unique(df$parsetting))
  {
    subdata2 <- df[df$parsetting == var,]; names(subdata2)
    laM <- unique(subdata2$sls.sim.lambda_M)
    muM <- unique(subdata2$sls.sim.mu_M)
    laS <- unique(subdata2$sls.sim.lambda_S)
    muS <- unique(subdata2$sls.sim.mu_S)
    conditioning <- unique(subdata2$sls.cond)
    pars_string <- paste0(laM,"-",muM,"-",laS,"-",muS,"-",conditioning)

    #lambda
    pdfname <- paste0("lambdaM_", pars_string)
    grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
    p5 <- ggplot2::ggplot(subdata2,
                          ggplot2::aes(x = sls.MLE.lambda_M,
                                       y = DDD.MLE.lambda_M
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
    ggplot2::labs(x = expression(lambda[M]^sls),
                  y = expression(lambda[M]^DDD)
    ) +
    ggplot2::theme_bw(
    ) +
    ggtitle(
      bquote(list(
        "Parameter setting: " ~
        lambda[M]==.(laM),
        mu[M]==.(muM),
        lambda[S]==.(laS),
        mu[S]==.(muS),
        cond==.(conditioning)
      ))
    ) +
    theme(axis.title = element_text(size = 18)); print(p5)
    grDevices::dev.off()

    #mu
    pdfname <- paste0("muM_", pars_string)
    grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
    p6 <- ggplot2::ggplot(subdata2,
                          ggplot2::aes(x = sls.MLE.mu_M,
                                       y = DDD.MLE.mu_M
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
    ggplot2::labs(x = expression(mu[M]^sls),
                  y = expression(mu[M]^DDD)
    ) +
    ggplot2::theme_bw(
    ) +
    ggtitle(
      bquote(list(
        "Parameter setting: " ~
        lambda[M]==.(laM),
        mu[M]==.(muM),
        lambda[S]==.(laS),
        mu[S]==.(muS),
        cond==.(conditioning)
      ))
    ) +
    theme(axis.title = element_text(size = 18)); print(p6)
    grDevices::dev.off()

    #lambda - mu
    pdfname <- paste0("lambdaM-muM_", pars_string)
    grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
    p7 <- ggplot2::ggplot(subdata2,
                          ggplot2::aes(x = (sls.MLE.lambda_M - sls.MLE.mu_M),
                                       y = (DDD.MLE.lambda_M - DDD.MLE.mu_M)
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
    ggplot2::labs(x = expression(lambda[M]^sls - mu[M]^sls),
                  y = expression(lambda[M]^DDD - mu[M]^DDD)
    ) +
    ggplot2::theme_bw(
    ) +
    ggtitle(
      bquote(list(
        "Parameter setting: " ~
        lambda[M]==.(laM),
        mu[M]==.(muM),
        lambda[S]==.(laS),
        mu[S]==.(muS),
        cond==.(conditioning)
      ))
    ) +
    theme(axis.title = element_text(size = 18)); print(p7)
    grDevices::dev.off()

    #mu/lambda
    #lambda - mu
    pdfname <- paste0("ratiomuMlambdaM_", pars_string)
    grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
    p8 <- ggplot2::ggplot(subdata2,
                          ggplot2::aes(x = (sls.MLE.mu_M/sls.MLE.lambda_M),
                                       y = (DDD.MLE.mu_M/DDD.MLE.lambda_M)
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
      ggplot2::labs(x =  expression(mu[M]^sls/lambda[M]^sls) ,
                    y =  expression(mu[M]^DDD/lambda[M]^DDD)
    ) +
    ggplot2::theme_bw(
    ) +
    ggtitle(
      bquote(list(
        "Parameter setting: " ~
        lambda[M]==.(laM),
        mu[M]==.(muM),
        lambda[S]==.(laS),
        mu[S]==.(muS),
        cond==.(conditioning)
      ))
    ) +
    theme(axis.title = element_text(size = 18)); print(p8)
    grDevices::dev.off()
  }

}

showplot.cloud <- function(data) {

  library(ggplot2); library(reshape2); library(scales)

  # home_dir <- substring(getwd(), 1, 21)
  # proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir, "/Progress/")))
  # folder_name <- paste0(home_dir, "/Progress/", list.files(paste0(home_dir, "/Progress"))[proj.coords])
  # results_folder <- paste0(folder_name, "/results")
  results_folder <- data$directory[1]
  path <- paste0(results_folder, "/cloudplots")

  if (!file.exists(path)) {dir.create(file.path(path), showWarnings = FALSE)}

  data2 <- data
  data2$parsetting <- interaction(data2$sim.lambda_M, data2$sim.mu_M,
                                  data2$sim.lambda_S, data2$sim.mu_S,
                                  data2$cond)

  ds <- data2[data2$model == "sls",]; names(ds) <- paste0("sls.", names(data2))
  dD <- data2[data2$model == "DDD",]; names(dD) <- paste0("DDD.", names(data2))
  df <- data.frame(ds, dD); names(df)
  df$parsetting <- interaction(df$sls.sim.lambda_M, df$sls.sim.mu_M,
                               df$sls.sim.lambda_S, df$sls.sim.mu_S,
                               df$sls.cond)

  for (var in unique(df$parsetting))
  {
    subdata2 <- df[df$parsetting == var,]; names(subdata2)
    laM <- unique(subdata2$sls.sim.lambda_M)
    muM <- unique(subdata2$sls.sim.mu_M)
    laS <- unique(subdata2$sls.sim.lambda_S)
    muS <- unique(subdata2$sls.sim.mu_S)
    conditioning <- unique(subdata2$sls.cond)
    pars_string  <- paste0(laM,"-",muM,"-",laS,"-",muS,"-",conditioning)

    #lambda vs mu
    pa <- ggplot2::ggplot(subdata2,
                          ggplot2::aes(x = sls.MLE.lambda_M,
                                       y = sls.MLE.mu_M
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
    ggplot2::labs(x = expression(lambda[M]^sls),
                  y = expression(mu[M]^sls)
    ) +
    ggplot2::theme_bw(
    ) +
    ggtitle(
      bquote(list(
        "Parameter setting: " ~
          lambda[M]==.(laM),
        mu[M]==.(muM),
        lambda[S]==.(laS),
        mu[S]==.(muS),
        cond==.(conditioning)
      ))
    ); #print(pa)

    #log lambda vs log mu
    pb <- ggplot2::ggplot(subdata2,
                          ggplot2::aes(x = sls.MLE.lambda_M,
                                       y = sls.MLE.mu_M
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
    ggplot2::labs(x = expression(lambda[M]^sls),
                  y = expression(mu[M]^sls)
    ) +
    ggplot2::theme_bw(
    ) +
    ggtitle(
      bquote(list(
        "Parameter setting: " ~
        lambda[M]==.(laM),
        mu[M]==.(muM),
        lambda[S]==.(laS),
        mu[S]==.(muS),
        cond==.(conditioning)
      ))
    ) +
    scale_x_continuous(trans = log2_trans()
    ) +
    scale_y_continuous(trans = log2_trans()
    ); #print(pb)

    #(lambda-mu) vs (mu/lambda)
    pc <- ggplot2::ggplot(subdata2,
                          ggplot2::aes(x = (sls.MLE.lambda_M - sls.MLE.mu_M),
                                       y = (sls.MLE.mu_M / sls.MLE.lambda_M)
                          )
    ) +
      ggplot2::geom_point(color = "firebrick4"
      ) +
      ggplot2::labs(x = expression(lambda[M]^sls - mu[M]^sls),
                    y = expression(mu[M]^sls / lambda[M]^sls)
      ) +
      ggplot2::theme_bw(
      ) +
      ggtitle(
        bquote(list(
          "Parameter setting: " ~
          lambda[M]==.(laM),
          mu[M]==.(muM),
          lambda[S]==.(laS),
          mu[S]==.(muS),
          cond==.(conditioning)
        ))
      ); #print(pc)

    #log(lambda-mu) vs log(mu/lambda)
    pd <- ggplot2::ggplot(subdata2,
                          ggplot2::aes(x = (sls.MLE.lambda_M - sls.MLE.mu_M),
                                       y = (sls.MLE.mu_M / sls.MLE.lambda_M)
                          )
    ) +
      ggplot2::geom_point(color = "firebrick4"
      ) +
      ggplot2::labs(x = expression(lambda[M]^sls - mu[M]^sls),
                    y = expression(mu[M]^sls / lambda[M]^sls)
      ) +
      ggplot2::theme_bw(
      ) +
      ggtitle(
        bquote(list(
          "Parameter setting: " ~
          lambda[M]==.(laM),
          mu[M]==.(muM),
          lambda[S]==.(laS),
          mu[S]==.(muS),
          cond==.(conditioning)
        ))
      ) +
      scale_x_continuous(trans = log2_trans()
      ) +
      scale_y_continuous(trans = log2_trans()
      ); #print(pd)

    #print outputs on pdf
    pdfname <- paste0("lambda_vs_mu", pars_string)
    grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
    print(pa)
    grDevices::dev.off()

    pdfname <- paste0("lambda_vs_mu[logscale]", pars_string)
    grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
    print(pb)
    grDevices::dev.off()

    pdfname <- paste0("lambda-mu_vs_muoverlambda", pars_string)
    grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
    print(pc)
    grDevices::dev.off()

    pdfname <- paste0("lambda-mu_vs_muoverlambda[logscale]", pars_string)
    grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
    print(pd)
    grDevices::dev.off()

  }

  }

# data analysis
data <- collect_data()
showplot.box(data)
showplot.correlation(data)
showplot.cloud(data)
