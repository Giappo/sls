#DATASETS TO TEST
data_folder <- paste0(getwd(), "/data/")
variable_name2string <- function(v1) {
  deparse(substitute(v1))
}

#pure branching: basic case
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6); mus <- c(0.3)
tb_when <- c(50); tb_who <- c(1); ts_when <- c(); ts_who <- c()
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_pure_branching1 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_pure_branching1, file = paste0(data_folder, variable_name2string(dataset_pure_branching1),".RData"))

#pure branching: second case
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6); mus <- c(0.3)
tb_when <- c(20, 40, 70, 90); tb_who <- c(1, 2, 3, 4); ts_when <- c(); ts_who <- c()
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/20; tf <- tf/20; ts[1,] <- ts[1,]/20; tb[1,] <- tb[1,]/20;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_pure_branching2 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_pure_branching2, file = paste0(data_folder, variable_name2string(dataset_pure_branching2),".RData"))

#pure branching: third case
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6); mus <- c(0.3)
tb_when <- c(20, 40, 45, 70, 80, 90); tb_who <- c(1, 2, 1, 3, 2, 4); ts_when <- c(); ts_who <- c()
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/20; tf <- tf/20; ts[1,] <- ts[1,]/20; tb[1,] <- tb[1,]/20;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_pure_branching3 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_pure_branching3, file = paste0(data_folder, variable_name2string(dataset_pure_branching3),".RData"))

#pure shifting: basic case
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4); mus <- c(0.3, 0.1)
tb_when <- c(); tb_who <- c(); ts_when <- c(50); ts_who <- c(1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_pure_shifting1 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_pure_shifting1, file = paste0(data_folder, variable_name2string(dataset_pure_shifting1),".RData"))

#pure shifting: second case
ti <- 0; tf <- 120; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4, 0.3, 0.2, 0.1); mus <- c(0.3, 0.2, 0.1, 0.05, 0.02)
tb_when <- c(); tb_who <- c(); ts_when <- c(20, 40, 70, 90); ts_who <- c(1,1,1,1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_pure_shifting2 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_pure_shifting2, file = paste0(data_folder, variable_name2string(dataset_pure_shifting2),".RData"))

#pure shifting: third case
ti <- 0; tf <- 300; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4, 0.3, 0.2, 0.1, 0.5, 0.3, 0.1, 0.05, 0.6); mus <- c(0.3, 0.2, 0.1, 0.05, 0.02, 0.3, 0.2, 0.1, 0.05, 0.02)
tb_when <- c(); tb_who <- c(); ts_when <- c(20, 40, 70, 90, 120, 150, 190, 210, 250); ts_who <- rep(1,length(ts_when))
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
den <- 50; ti <- ti/den; tf <- tf/den; ts[1,] <- ts[1,]/den; tb[1,] <- tb[1,]/den;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_pure_shifting3 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_pure_shifting3, file = paste0(data_folder, variable_name2string(dataset_pure_shifting3),".RData"))

#Fork's example
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4); mus <- c(0.3, 0.1)
tb_when <- c(60); tb_who <- c(1); ts_when <- c(30); ts_who <- c(1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_Fork <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_Fork, file = paste0(data_folder, variable_name2string(dataset_Fork),".RData"))

#R's example
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4); mus <- c(0.3, 0.1)
tb_when <- c(30); tb_who <- c(1); ts_when <- c(60); ts_who <- c(1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_Rampal <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_Rampal, file = paste0(data_folder, variable_name2string(dataset_Rampal),".RData"))

#B's example
ti <- 0; tf <- 150; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4, 0.2); mus <- c(0.3, 0.1, 0.05)
tb_when <- c(60); tb_who <- c(1); ts_when <- c(30,100); ts_who <- c(1,1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_Bart <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_Bart, file = paste0(data_folder, variable_name2string(dataset_Bart),".RData"))

# #B's example (shorter times)
# ti <- 0; tf <- 150; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
# lambdas <- c(0.6, 0.4, 0.2); mus <- c(0.3, 0.1, 0.05)
# tb_when <- c(60); tb_who <- c(1); ts_when <- c(30,100); ts_who <- c(1,1)
# tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
# common_factor <- 20; ti <- ti/common_factor; tf <- tf/common_factor; ts[1,] <- ts[1,]/common_factor; tb[1,] <- tb[1,]/common_factor;
# times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
# dataset_Bart2 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
# save(dataset_Bart2, file = paste0(data_folder, variable_name2string(dataset_Bart2),".RData"))

#one branching and three shifts, one before it and two after
ti <- 0; tf <- 15; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.5, 0.4, 0.2, 0.35); mus <- c(0.2, 0.1, 0.05, 0.15)
tb_when <- c(3); tb_who <- c(1); ts_when <- c(1, 6, 10); ts_who <- c(1, 1, 2)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_1 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_1, file = paste0(data_folder, variable_name2string(dataset_1),".RData"))

#most complicated example so far: it didn't work (yet) :(
ti <- 0; tf <- 21; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.5, 0.4, 0.2, 0.35, 0.25); mus <- c(0.2, 0.1, 0.05, 0.15, 0.05)
tb_when <- c(3, 10, 12, 16); tb_who <- c(1, 1, 3, 3); ts_when <- c(1, 6, 15, 19); ts_who <- c(1, 2, 3, 1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_2 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_2, file = paste0(data_folder, variable_name2string(dataset_2),".RData"))

#most complicated example so far: just bigger time points
ti <- 0; tf <- 21; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.5, 0.4, 0.2, 0.35, 0.25); mus <- c(0.2, 0.1, 0.05, 0.15, 0.05)
tb_when <- c(3, 10, 12, 16); tb_who <- c(1, 1, 3, 3); ts_when <- c(1, 6, 15, 19); ts_who <- c(1, 2, 3, 1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/5; tf <- tf/5; ts[1,] <- ts[1,]/5; tb[1,] <- tb[1,]/5;
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_3 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_3, file = paste0(data_folder, variable_name2string(dataset_3),".RData"))

# #most complicated example so far: just even bigger time points
# ti <- 0; tf <- 21; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
# lambdas <- c(0.5, 0.4, 0.2, 0.35, 0.25); mus <- c(0.2, 0.1, 0.05, 0.15, 0.05)
# tb_when <- c(3, 10, 12, 16); tb_who <- c(1, 1, 3, 3); ts_when <- c(1, 6, 15, 19); ts_who <- c(1, 2, 3, 1)
# tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
#
# times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
# dataset_4 <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
# save(dataset_4, file = paste0(data_folder, variable_name2string(dataset_4),".RData"))

#dataset to test pure birth: branching, shift, branching, branching
ti <- 0; tf <- 5; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.2, 0.1); mus <- rep(0,length(lambdas))
tb_when <- c(1, 3, 4); tb_who <- c(1, 1, 1); ts_when <- c(2); ts_who <- c(2)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
times_matrix <- arrange_times_matrix(ti = ti, tb = tb, ts = ts, tf = tf)
dataset_PB <- list(times_matrix = times_matrix, lambdas = lambdas, mus = mus)
save(dataset_PB, file = paste0(data_folder, variable_name2string(dataset_PB),".RData"))

#clean
rm(ti,ts,tb,tf,lambdas,mus, tb_when, tb_who, ts_when, ts_who, N0, den, variable_name2string, data_folder, times_matrix, common_factor)
