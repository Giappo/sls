maxsims <- 50
idparsopt <- c(1:2)
pars <- c(0.3, 0.1, 0.6, 0.05)
t0 <- c(10, 4)

Npars <- length(pars)
idparsfix <- (1:Npars)[-idparsopt]
initparsopt <- c(0.2, 0.1, 0.2, 0.1)[idparsopt]
parsfix <- pars[idparsfix]
MLE <- matrix(NA, nrow = maxsims, ncol = (Npars + 3))
cond <- 3
sims <- vector("list", maxsims)
for (s in 1:maxsims) #19 is critical
{
  set.seed(s)

  sims[[s]] <- sls::sls_sim(pars = pars, t0 = t0, cond = cond)
  brtsM <- sims[[s]]$brts[[1]]; brtsS <- sims[[s]]$brts[[2]]; brtsM; brtsS

  # DDD_test <- sls::DDD_sim(pars = pars, t0 = t0, starting_species = c(2, 1), cond = 1)
  # L <- DDD_test$L
  # L0 <- L[L[,5] == 0, 1:5]; dim(L0) <- c(sum(L[,5] == 0), 5); L0
  # brtsM <- c(t0[1], DDD:::L2brts(L0))
  # L1 <- L[L[,5] == 1, 1:5]; dim(L1) <- c(sum(L[,5] == 1), 5); L1
  # if(nrow(L1) > 1) {brtsS <- c(t0[2], DDD:::L2brts(L1))}else{brtsS <- c(t0[2], L1[1,1])}; print(brtsS)

  MLE[s,] <- unlist(
    sls::sls_ML(initparsopt = initparsopt, parsfix = parsfix, idparsopt = idparsopt, idparsfix = idparsfix,
                cond = cond, brtsM = brtsM, brtsS = brtsS)
  )
}; MLE
fails = which(MLE[,1] == -1)
how_many_fails <- length(fails); how_many_fails
MLE_clear <- MLE[MLE[,1] != -1,]
apply(MLE_clear[,idparsopt], MARGIN = 2, FUN = function(x) median(x[x != -1])); pars[idparsopt]
df <- as.data.frame(MLE_clear); names(df) <- c("lambda_M", "mu_M", "lambda_S", "mu_S", "LL", "df", "conv"); head(df)
histograms <- list(length(idparsopt))
for (i in seq_along(idparsopt))
{
  histograms[[i]] <- ggplot2::qplot(df[,idparsopt[i]], geom = 'histogram') +
    ggplot2::geom_vline(xintercept = pars[idparsopt[i]], col = "red") +
    ggplot2::ggtitle(paste0("Results for ", names(df)[idparsopt[i]], " over ", maxsims - how_many_fails, "/", maxsims, " trees")) +
    ggplot2::xlab(names(df)[idparsopt[i]])
}
histograms[[1]]
histograms[[2]]

#problematic trees
i <- 1
i <- i + 1; sims[[fails[i]]]

L1 <- sims[[fails[i]]]$L[[1]]
DDD::L2brts(L1)
plot(DDD::L2phylo(L1))
