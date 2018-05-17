#DATASETS TO TEST

#pure branching
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6); mus <- c(0.3)
tb_when <- c(50); tb_who <- c(1); ts_when <- c(); ts_who <- c()
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_pure_branching <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#pure shifting
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4); mus <- c(0.3, 0.1)
tb_when <- c(); tb_who <- c(); ts_when <- c(50); ts_who <- c(1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_pure_shifting <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#R's example
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4); mus <- c(0.3, 0.1)
tb_when <- c(30); tb_who <- c(1); ts_when <- c(60); ts_who <- c(1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_Rampal <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#B's example
ti <- 0; tf <- 150; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4, 0.2); mus <- c(0.3, 0.1, 0.05)
tb_when <- c(60); tb_who <- c(1); ts_when <- c(30,100); ts_who <- c(1,1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_Bart <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#one branching and three shifts, one before it and two after
ti <- 0; tf <- 15; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.5, 0.4, 0.2, 0.35); mus <- c(0.2, 0.1, 0.05, 0.15)
tb_when <- c(3); tb_who <- c(1); ts_when <- c(1, 6, 10); ts_who <- c(1, 1, 2)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_1 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#most complicated example so far: it didn't work (yet) :(
ti <- 0; tf <- 21; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.5, 0.4, 0.2, 0.35, 0.25); mus <- c(0.2, 0.1, 0.05, 0.15, 0.05)
tb_when <- c(3, 10, 12, 16); tb_who <- c(1, 1, 3, 3); ts_when <- c(1, 6, 15, 19); ts_who <- c(1, 2, 3, 1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_2 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#most complicated example so far: just bigger time points
ti <- 0; tf <- 21; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.5, 0.4, 0.2, 0.35, 0.25); mus <- c(0.2, 0.1, 0.05, 0.15, 0.05)
tb_when <- c(3, 10, 12, 16); tb_who <- c(1, 1, 3, 3); ts_when <- c(1, 6, 15, 19); ts_who <- c(1, 2, 3, 1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/5; tf <- tf/5; ts[1,] <- ts[1,]/5; tb[1,] <- tb[1,]/5;
dataset_3 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#most complicated example so far: just even bigger time points
ti <- 0; tf <- 21; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.5, 0.4, 0.2, 0.35, 0.25); mus <- c(0.2, 0.1, 0.05, 0.15, 0.05)
tb_when <- c(3, 10, 12, 16); tb_who <- c(1, 1, 3, 3); ts_when <- c(1, 6, 15, 19); ts_who <- c(1, 2, 3, 1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);

dataset_4 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#clean
rm(ti,ts,tb,tf,lambdas,mus, tb_when, tb_who, ts_when, ts_who, N0)
