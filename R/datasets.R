# RJCB: suggest to put this in a function, e.g. create_dataset
#DATASETS TO TEST

#pure branching: basic case
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6); mus <- c(0.3)
tb_when <- c(50); tb_who <- c(1); ts_when <- c(); ts_who <- c()
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_pure_branching1 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#pure branching: second case
ti <- 0; tf <- 120; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6); mus <- c(0.3)
tb_when <- c(20, 40, 70, 90); tb_who <- c(1,2,3,4); ts_when <- c(); ts_who <- c()
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_pure_branching2 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#pure branching: third case
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6); mus <- c(0.3)
tb_when <- c(20, 40, 70); tb_who <- c(1, 2, 3); ts_when <- c(); ts_who <- c()
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/20; tf <- tf/20; ts[1,] <- ts[1,]/20; tb[1,] <- tb[1,]/20;
dataset_pure_branching3 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#pure shifting: basic case
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4); mus <- c(0.3, 0.1)
tb_when <- c(); tb_who <- c(); ts_when <- c(50); ts_who <- c(1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_pure_shifting1 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#pure shifting: advanced case
ti <- 0; tf <- 120; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4, 0.3, 0.2, 0.1); mus <- c(0.3, 0.2, 0.1, 0.05, 0.02)
tb_when <- c(); tb_who <- c(); ts_when <- c(20, 40, 70, 90); ts_who <- c(1,1,1,1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_pure_shifting2 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#pure shifting: advanced case 2
ti <- 0; tf <- 300; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4, 0.3, 0.2, 0.1, 0.5, 0.3, 0.1, 0.05, 0.6); mus <- c(0.3, 0.2, 0.1, 0.05, 0.02, 0.3, 0.2, 0.1, 0.05, 0.02)
tb_when <- c(); tb_who <- c(); ts_when <- c(20, 40, 70, 90, 120, 150, 190, 210, 250); ts_who <- rep(1,length(ts_when))
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
den <- 50; ti <- ti/den; tf <- tf/den; ts[1,] <- ts[1,]/den; tb[1,] <- tb[1,]/den;
dataset_pure_shifting3 <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#Fork's example
ti <- 0; tf <- 100; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.6, 0.4); mus <- c(0.3, 0.1)
tb_when <- c(60); tb_who <- c(1); ts_when <- c(30); ts_who <- c(1)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
ti <- ti/10; tf <- tf/10; ts[1,] <- ts[1,]/10; tb[1,] <- tb[1,]/10;
dataset_Fork <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

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

#dataset to test pure birth: branching, shift, branching, branching
ti <- 0; tf <- 5; N0 <- 1; ts <- tb <- matrix(NA, nrow = 2);
lambdas <- c(0.2, 0.1); mus <- rep(0,length(lambdas))
tb_when <- c(1, 3, 4); tb_who <- c(1, 1, 1); ts_when <- c(2); ts_who <- c(2)
tb <- rbind(tb_when, tb_who); ts <- rbind(ts_when, ts_who);
dataset_PB <- list(ti = ti, tb = tb, ts = ts, tf = tf, lambdas = lambdas, mus = mus)

#clean
rm(ti,ts,tb,tf,lambdas,mus, tb_when, tb_who, ts_when, ts_who, N0, den)
