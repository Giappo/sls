# rm(list = ls())
# FUNCTIONS -----------------------------------------------------------------------------
# functions - utility -----------------------------------------------------------------------------

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

correct_model_label <- function(old_model_name = "sls", new_model_name = "slsP") {

  db_dir <- get.dropbox.folder()
  home_dir <- db_dir
  home_dir <- paste0(db_dir, "\\university\\Progress\\"); list.files(home_dir)
  proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir)))
  folder_name <- paste0(home_dir, list.files(paste0(home_dir))[proj.coords])
  results_mother_folder <- paste0(folder_name, "/results")
  # results_folder <- choose.dir(default = results_mother_folder)
  results_folder <- results_mother_folder
  datasets <- list.files(results_folder, pattern = "^[0]")
  Nd <- length(datasets)

  d <- 1
  for (d in 1:Nd)
  {# d loop
    print(d)

    local_path <- paste0(results_folder, "//", datasets[d]);
    setwd(local_path)
    names0 <- list.files(pattern = "[.]txt", path = local_path, full.names = TRUE)
    names1 <- gsub(".*/","", names0)
    model_names <- unique(gsub("_MLE.*", "", names1))
    if (old_model_name %in% model_names)
    {
      filenames <- paste0(old_model_name, "_MLE")
      files <- list.files(pattern = paste0(filenames), path = local_path, full.names = FALSE); length(files)
      s <- 1
      for (s in 1:length(files))
      {# s loop
        new_name <- gsub(old_model_name, new_model_name, files[s])
        file.rename(from = files[s], to = new_name)
      }# s loop
    }
  }# d loop
}

delete_useless_files <- function() {

  db_dir <- get.dropbox.folder()
  home_dir <- db_dir
  home_dir <- paste0(db_dir, "\\university\\Progress\\"); list.files(home_dir)
  proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir)))
  folder_name <- paste0(home_dir, list.files(paste0(home_dir))[proj.coords])
  results_mother_folder <- paste0(folder_name, "/results")
  # results_folder <- choose.dir(default = results_mother_folder)
  results_folder <- results_mother_folder
  datasets <- list.files(results_folder, pattern = "^[0]")
  Nd <- length(datasets)

  d <- 1
  for (d in 1:Nd)
  {# d loop
    print(d)

    local_path <- paste0(results_folder, "//", datasets[d]);
    setwd(local_path)
    Sys.chmod(local_path, "777")
    files_MLE   <- list.files(pattern = "[.]txt", path = local_path, full.names = TRUE)
    files_slurm <- list.files(pattern = "slurm", path = local_path, full.names = TRUE)
    files_all   <- list.files(path = local_path, full.names = TRUE)
    not_in <- function(big_set, sub_set) {big_set[!('%in%'(big_set, sub_set))]}
    files1 <- not_in(big_set = files_all, sub_set = files_MLE)
    files2 <- not_in(big_set = files1   , sub_set = files_slurm)

    names1 <- gsub(".*/","", files2)
    if (length(files2) > 0)
    {
      s <- 1
      for (s in 1:length(files2))
      {# s loop
        file.remove(files2[s])
      }# s loop
    }
  }# d loop

}

collect_data <- function() {

  db_dir <- get.dropbox.folder()
  home_dir <- db_dir
  home_dir <- paste0(db_dir, "\\university\\Progress\\"); list.files(home_dir)
  proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir)))
  folder_name <- paste0(home_dir, list.files(paste0(home_dir))[proj.coords])
  results_mother_folder <- paste0(folder_name, "/results")
  # results_folder <- choose.dir(default = results_mother_folder)
  results_folder <- results_mother_folder
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
    model_names <- unique(gsub("_MLE.*", "", names1))
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
        # right_ids <- results2[results2[, 1] != -1, "tree_id"]; length(right_ids)
        # results3  <- results2[results2[, "tree_id"] %in% right_ids,]; dim(results3)
        results3  <- results2
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
      # model_ids[[m]] <- right_ids
      d_results <- rbind(d_results, results5)
      # rm(results0, results1, results2, results3, results4, results5,
      #    right_ids, files, filenames)
    }# m loop
    # d_ids <- Reduce(intersect, model_ids); d_ids <- sort(d_ids)
    # d_results <- d_results[(d_results[,"tree_id"] %in% d_ids),]
    all_results <- rbind(all_results, d_results)
    rm(d_results)

  }# d loop

  # for (m in seq_along(model_names)[-1])
  # {
  #   testit::assert(
  #   dim(all_results[all_results$model == model_names[m - 1],]) == dim(all_results[all_results$model == model_names[m],])
  #   )
  # }

  all_results$directory <- results_folder

  return(all_results)
}

data_overview <- function(data) {

  data1 <- data
  all_models <- levels(droplevels(unique(data1$model)))
  cat(paste0("There are results for the following models: ")); print(all_models)
  simpars.names <- c("sim.lambda_M","sim.mu_M","sim.lambda_S","sim.mu_S","cond")
  data2 <- data1[, simpars.names]
  parsettings <- unique(data2); Nd <- nrow(parsettings); Nd
  Ndm <- matrix(0, ncol = length(all_models), nrow = Nd)
  d_models <- vector("list", Nd); d <- 1
  for (d in 1:Nd)
  {# d loop
    primed_data <- 2^data2[,1]*3^data2[,2]*5^data2[,3]*7^data2[,4]*11^data2[,5]
    primed_pars <- 2^parsettings[d,1]*3^parsettings[d,2]*5^parsettings[d,3]*7^parsettings[d,4]*11^parsettings[d,5]
    drows <- which(primed_data == primed_pars)
    data3 <- data1[drows,]
    d_models[[d]] <- levels(droplevels(unique(data3$model)))
    m <- 1
    for (m in seq_along(all_models))
    {
      Ndm[d, m] <- sum(data3$model == all_models[m])
    }
  }; d_models
  overview <- matrix(NA, nrow = Nd, ncol = (ncol(parsettings) + length(all_models)))
  colnames(overview) <- c(colnames(data2), all_models)
  # overview[,1:ncol(parsettings)] <- parsettings
  overview[, simpars.names] <- as.matrix(parsettings)

  # for (d in 1:Nd)
  # {# d loop
  #   overview[d, (length(simpars.names) + 1): ncol(overview)] <- d_models[[d]] %in% all_models
  # }
  overview[, (length(simpars.names) + 1): ncol(overview)] <- Ndm
  return(overview)
}

load_trees <- function(s) {

  testit::assert(is.numeric(s))
  testit::assert(floor(s) == ceiling(s))

  db_dir <- get.dropbox.folder()
  home_dir <- db_dir
  home_dir <- paste0(db_dir, "\\university\\Progress\\"); list.files(home_dir)
  proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir)))
  folder_name <- paste0(home_dir, list.files(paste0(home_dir))[proj.coords])
  results_mother_folder <- paste0(folder_name, "/results")
  # results_folder <- choose.dir(default = results_mother_folder)
  results_folder <- results_mother_folder
  datasets <- list.files(results_folder, pattern = "^[0]")

  dataset_pars <- vector("list", Nd <- length(datasets))
  for (d in 1:Nd)
  {
    dataset_pars[[d]] <- as.numeric(unlist( strsplit(x = datasets[d], split = "-")))
  }
  cat("Available datasets are: ")
  print(datasets)
  x <- readline("What dataset do you choose? \n")

  dataset_folder <- paste0(results_folder, "/" ,datasets[as.numeric(x)])
  data_folder <- paste0(dataset_folder, "/", "data")
  # list.files(data_folder)
  if (missing(s)) {
    s <- readline("What tree do you choose? \n")
  }
  tree_file <- paste0(data_folder, "/", "sim_", s ,".Rdata")
  # simtrees <- NULL
  if (length(s) == 1)
    {
    while (!file.exists(tree_file)) {
      print("This tree is not present. \n")
      s <- readline("What tree do you choose? \n")
      tree_file <- paste0(data_folder, "/", "sim_", s ,".Rdata")
    }
    simtrees <- get(load(tree_file))
  }else
  {
    simtrees <- list(); ss <- 1
    for (ss in s)
    {
      if (!file.exists(tree_file[ss]))
      {
        simtrees[[ss]] <- NULL
      }else
      {
        simtrees[[ss]] <- get(load(tree_file[ss]))
      }
    }
  }

  return(simtrees)
}

allow_model_comparison <- function(data, model1, model2) {
  models <- unique(data$model)
  if (!(model1 %in% models)) {stop('model1 is not present in the data')}
  if (!(model2 %in% models)) {stop('model2 is not present in the data')}
  data1 <- data
  models <- c(model1, model2)

  data2 <- data1[,c("sim.lambda_M","sim.mu_M","sim.lambda_S","sim.mu_S","cond")]
  parsettings <- unique(data2); Nd <- nrow(parsettings); Nd
  all_results <- NULL; d <- 1
  for (d in 1:Nd)
  {# d loop
    primed_data <- 2^data2[,1]*3^data2[,2]*5^data2[,3]*7^data2[,4]*11^data2[,5]
    primed_pars <- 2^parsettings[d,1]*3^parsettings[d,2]*5^parsettings[d,3]*7^parsettings[d,4]*11^parsettings[d,5]
    drows <- which(primed_data == primed_pars)
    data3 <- data1[drows,]

    d_results  <- NULL; model_ids <- vector("list", length(models)); m <- 1
    for (m in seq_along(models))
    {# m loop
      model     <- models[m]
      mrows     <- data3$model == model
      if (sum(mrows) != 0)
      {
        results2  <- data3[mrows,]
        right_ids <- results2[results2[, 1] != -1, "tree_id"]; length(right_ids)
        results3  <- results2[results2[, "tree_id"] %in% right_ids,]; dim(results3)
        results4  <- results3[order(results3$tree_id),]

        testit::assert(results4$cond == floor(results4$cond) | results4$cond == 0)
        model_ids[[m]] <- right_ids
        d_results <- rbind(d_results, results4)

        rm(results2, results3, results4)
      }
    }# m loop
    d_ids <- Reduce(intersect, model_ids); d_ids <- sort(d_ids)
    d_results <- d_results[(d_results[,"tree_id"] %in% d_ids),]
    all_results <- rbind(all_results, d_results)
    rm(d_results)
  }# d loop

  testit::assert(
    sum(all_results$model == model1) == sum(all_results$model == model2)
  )

  return(all_results)
}

# functions - data analysis -----------------------------------------------------------------------------

showplot.box <- function(data) {

  library(ggplot2); library(reshape2); library(scales)

  # home_dir <- substring(getwd(), 1, 21)
  # proj.coords <- pmatch(x = "RQ4", table = list.files(paste0(home_dir, "/Progress/")))
  # folder_name <- paste0(home_dir, "/Progress/", list.files(paste0(home_dir, "/Progress"))[proj.coords])
  # results_folder <- paste0(folder_name, "/results")
  results_folder <- data$directory[1]
  models <- unique(data$model)
  models <- sort(models)
  testit::assert(length(models) == 2)
  model_comparison_folder <- paste0(results_folder, "\\", models[1], "_vs_", models[2])
  if (!file.exists(model_comparison_folder)) {dir.create(model_comparison_folder, showWarnings = FALSE)}
  path <- paste0(model_comparison_folder, "/boxplots")
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
  models <- unique(data$model)
  models <- sort(models)
  testit::assert(length(models) == 2)
  model_comparison_folder <- paste0(results_folder, "\\", models[1], "_vs_", models[2])
  if (!file.exists(model_comparison_folder)) {dir.create(model_comparison_folder)}
  path <- paste0(model_comparison_folder, "/correlationplots")
  if (!file.exists(path)) {dir.create(file.path(path), showWarnings = FALSE)}

  data2 <- data; names(data2)
  data2$parsetting <- interaction(data2$sim.lambda_M, data2$sim.mu_M,
                                  data2$sim.lambda_S, data2$sim.mu_S,
                                  data2$cond)

  # ds <- data2[data2$model == models[1],]; names(ds) <- paste0(models[1], ".", names(data2)); names(ds)
  # dD <- data2[data2$model == models[2],]; names(dD) <- paste0(models[2], ".", names(data2)); names(dD)
  ds <- data2[data2$model == models[1],]; names(ds) <- paste0("mod1.", names(data2)); names(ds)
  dD <- data2[data2$model == models[2],]; names(dD) <- paste0("mod2.", names(data2)); names(dD)
  df <- data.frame(ds, dD); names(df)

  # sim.names <- names(df)[
  #   grepl(
  #   paste0(models[2],".sim"),
  #   names(df)
  #   ) |
  #   grepl(
  #     paste0(models[2],".cond"),
  #     names(df)
  #   )
  # ]
  # model1.names <- names(df)[
  #   grepl(
  #     paste0(models[1],".MLE"),
  #     names(df)
  #   )
  # ]
  # model2.names <- names(df)[
  #   grepl(
  #     paste0(models[2],".MLE"),
  #     names(df)
  #   )
  # ]
  # df$parsetting <- interaction(df[sim.names[1:length(sim.names)]])
  df$parsetting <- interaction(df$mod1.sim.lambda_M, df$mod1.sim.mu_M,
                               df$mod1.sim.lambda_S, df$mod1.sim.mu_S,
                               df$mod1.cond)

  for (var in unique(df$parsetting))
  {
    subdata2 <- df[df$parsetting == var,]; names(subdata2)
    laM <- unique(subdata2$mod1.sim.lambda_M)
    muM <- unique(subdata2$mod1.sim.mu_M)
    laS <- unique(subdata2$mod1.sim.lambda_S)
    muS <- unique(subdata2$mod1.sim.mu_S)
    conditioning <- unique(subdata2$mod1.cond)

    # laM <- unique(subdata2[sim.names[1]])
    # muM <- unique(subdata2[sim.names[2]])
    # laS <- unique(subdata2[sim.names[3]])
    # muS <- unique(subdata2[sim.names[4]])
    # conditioning <- unique(subdata2[sim.names[5]])

    pars_string <- paste0(laM,"-",muM,"-",laS,"-",muS,"-",conditioning)
    labA <- levels(droplevels(models[1]))
    labB <- levels(droplevels(models[2]))

    #lambda
    pdfname <- paste0("lambdaM_", pars_string)
    grDevices::pdf(file = paste0(path, "/", pdfname, ".pdf"));
    p5 <- ggplot2::ggplot(subdata2,
                          ggplot2::aes(x = mod1.MLE.lambda_M,
                                       y = mod2.MLE.lambda_M
                          # ggplot2::aes(x = eval(model1.names[1]),
                          #              y = eval(model2.names[1])
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
    ggplot2::labs(
      # x = expression(lambda[M]^models[1]),
      # y = expression(lambda[M]^models[2])
      x = bquote(lambda[M]^.(labA)),
      y = bquote(lambda[M]^.(labB))
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
                          ggplot2::aes(
                            x = mod1.MLE.mu_M,
                            y = mod2.MLE.mu_M
                          # ggplot2::aes(x = sls.MLE.mu_M,
                          #              y = DDD.MLE.mu_M
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
    ggplot2::labs(
      # x = expression(mu[M]^sls),
      # y = expression(mu[M]^DDD)
      x = bquote(mu[M]^.(labA)),
      y = bquote(mu[M]^.(labB))
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
                          ggplot2::aes(
                            # x = (sls.MLE.lambda_M - sls.MLE.mu_M),
                            # y = (DDD.MLE.lambda_M - DDD.MLE.mu_M)
                            x = mod1.MLE.lambda_M - mod1.MLE.mu_M,
                            y = mod2.MLE.lambda_M - mod2.MLE.mu_M
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
    ggplot2::labs(
      # x = expression(lambda[M]^sls - mu[M]^sls),
      # y = expression(lambda[M]^DDD - mu[M]^DDD)
      x = bquote(lambda[M]^.(labA) - mu[M]^.(labA)),
      y = bquote(lambda[M]^.(labB) - mu[M]^.(labB))
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
                          ggplot2::aes(
                            # x = (sls.MLE.mu_M/sls.MLE.lambda_M),
                            # y = (DDD.MLE.mu_M/DDD.MLE.lambda_M)
                            x = mod1.MLE.mu_M/mod1.MLE.lambda_M,
                            y = mod2.MLE.mu_M/mod2.MLE.lambda_M
                          )
    ) +
    ggplot2::geom_point(color = "firebrick4"
    ) +
      ggplot2::labs(
        # x =  expression(mu[M]^sls/lambda[M]^sls) ,
        # y =  expression(mu[M]^DDD/lambda[M]^DDD)
        x = bquote(mu[M]^.(labA)/lambda[M]^.(labA)),
        y = bquote(mu[M]^.(labB)/lambda[M]^.(labB))
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

  results_folder <- data$directory[1]
  models <- unique(data$model)
  models <- sort(models)

  m <- 1
  for (m in seq_along(models))
  { # m loop
    model_folder <- paste0(results_folder, "\\", models[m])
    if (!file.exists(model_folder)) {dir.create(model_folder)}
    path <- paste0(model_folder, "/cloudplots")
    if (!file.exists(path)) {dir.create(file.path(path), showWarnings = FALSE)}

    data2 <- data[data$model == models[m],]; names(data2)
    data2$parsetting <- interaction(data2$sim.lambda_M, data2$sim.mu_M,
                                    data2$sim.lambda_S, data2$sim.mu_S,
                                    data2$cond)

    df <- data2
    df$parsetting <- interaction(df$sim.lambda_M, df$sim.mu_M,
                                 df$sim.lambda_S, df$sim.mu_S,
                                 df$cond)

    for (var in unique(df$parsetting))
    { # parsetting loop
      subdata2 <- df[df$parsetting == var,]; names(subdata2)
      laM <- unique(subdata2$sim.lambda_M)
      muM <- unique(subdata2$sim.mu_M)
      laS <- unique(subdata2$sim.lambda_S)
      muS <- unique(subdata2$sim.mu_S)
      conditioning <- unique(subdata2$cond)

      pars_string <- paste0(laM,"-",muM,"-",laS,"-",muS,"-",conditioning)
      labA <- levels(droplevels(models[m]))

      #lambda vs mu
      pa <- ggplot2::ggplot(subdata2,
                            ggplot2::aes(
                              x = MLE.lambda_M,
                              y = MLE.mu_M
                            )
      ) +
      ggplot2::geom_point(color = "firebrick4"
      ) +
      ggplot2::labs(
        # x = expression(lambda[M]^sls),
        # y = expression(mu[M]^sls)
        x = bquote(lambda[M]^.(labA)),
        y = bquote(mu[M]^.(labA))
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
      ); pa#print(pa)

      #log lambda vs log mu
      pb <- ggplot2::ggplot(subdata2,
                            ggplot2::aes(
                              x = MLE.lambda_M,
                              y = MLE.mu_M
                            )
      ) +
      ggplot2::geom_point(color = "firebrick4"
      ) +
      ggplot2::labs(
        # x = expression(lambda[M]^sls),
        # y = expression(mu[M]^sls)
        x = bquote(lambda[M]^.(labA)),
        y = bquote(mu[M]^.(labA))
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
      ); pb#print(pb)

      #(lambda-mu) vs (mu/lambda)
      pc <- ggplot2::ggplot(subdata2,
                            ggplot2::aes(
                              x = (MLE.lambda_M - MLE.mu_M),
                              y = (MLE.mu_M / MLE.lambda_M)
                            )
      ) +
        ggplot2::geom_point(color = "firebrick4"
        ) +
        ggplot2::labs(
          # x = expression(lambda[M]^sls - mu[M]^sls),
          # y = expression(mu[M]^sls / lambda[M]^sls)
          x = bquote(lambda[M]^.(labA) - mu[M]^.(labA)),
          y = bquote(mu[M]^.(labA) / lambda[M]^.(labA))
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
        ); pc#print(pc)

      #log(lambda-mu) vs log(mu/lambda)
      pd <- ggplot2::ggplot(subdata2,
                            ggplot2::aes(
                              # x = (sls.MLE.lambda_M - sls.MLE.mu_M),
                              # y = (sls.MLE.mu_M / sls.MLE.lambda_M)
                              x = (MLE.lambda_M - MLE.mu_M),
                              y = (MLE.mu_M / MLE.lambda_M)
                            )
      ) +
        ggplot2::geom_point(color = "firebrick4"
        ) +
        ggplot2::labs(
          # x = expression(lambda[M]^sls - mu[M]^sls),
          # y = expression(mu[M]^sls / lambda[M]^sls)
          x = bquote(lambda[M]^.(labA) - mu[M]^.(labA)),
          y = bquote(mu[M]^.(labA) / lambda[M]^.(labA))
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
        ); pd#print(pd)

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
    } # parsetting loop
  }# m loop

}

# DATA ANALYSIS -----------------------------------------------------------------------------

data0 <- collect_data()
data_overview(data0)

# correct_model_label(old_model_name = "sls", new_model_name = "slsP")
# delete_useless_files()
# load_trees(1:50)

data01  <- allow_model_comparison(data0, model1 = "slsP", model2 = "DDD")
showplot.box(data01)
showplot.correlation(data01)
showplot.cloud(data01)



test <- load_trees(1:20)


