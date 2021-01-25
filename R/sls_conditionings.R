#' @title sls-conditioning
#' @author Giovanni Laudanno
#' @description Calculates three different kind of conditioning probabilities
#' @inheritParams default_params_doc
#' @return Three different conditioning probabilities
#' @export
pc_1shift <- function(
  pars_m,
  pars_s,
  brts_m,
  brts_s,
  cond,
  n_max = 1e2,
  n_0 = 2
) {
  lambdas <- c(pars_m[1], pars_s[1])
  mus     <- c(pars_m[2], pars_s[2])

  brts_m1 <- sort(abs(brts_m), decreasing = TRUE)
  brts_s1 <- sort(abs(brts_s), decreasing = TRUE)
  t_s <- brts_s1[1]
  t_c <- brts_m1[1]
  t_p <- 0
  aa <- abs(t_s - t_c); bb <- abs(t_p - t_s)

  if (n_0 != 2) {
    stop("Pc can be calculated only if phylogeny starts with a crown!")
  }

  p_s <- sls::p_t(
    t = bb,
    lambda = lambdas[2],
    mu = mus[2]
  )

  # the shift happens on the right (r) branch
  # legend: "r" is right; "l" is left
  # legend: "c" is crown; "s" is shift; "p" is present
  # legend: "aa" is t_s - t_c; "bb" is t_p - t_s
  # legend: "ns" is the number of species at the shift point
  # legend: "n_l" is "ns" for left branch; "n_r" is "ns" for right branch;

  nvec <- 1:n_max
  n_r <- row(matrix(NA, nrow = n_max, ncol = n_max))
  n_l <- col(matrix(NA, nrow = n_max, ncol = n_max))
  p_a_1 <- sls::p_t(t = aa, lambda = lambdas[1], mu = mus[1]); p_a_1
  u_a_1 <- sls::ut(t = aa, lambda = lambdas[1], mu = mus[1]); u_a_1
  p_b_1 <- sls::p_t(t = bb, lambda = lambdas[1], mu = mus[1]); p_b_1
  one_minus_p_b_1  <-
    sls::one_minus_pt(t = bb, lambda = lambdas[1], mu = mus[1])
  p_b_2 <- sls::p_t(t = bb, lambda = lambdas[2], mu = mus[2]); p_b_2
  p_n_r <- sls::pn(n = n_r, t = aa, lambda = lambdas[1], mu = mus[1])
  rownames(p_n_r) <- paste0("n_r=", nvec)
  colnames(p_n_r) <- paste0("n_l=", nvec)
  p_n_l <- sls::pn(n = n_l, t = aa, lambda = lambdas[1], mu = mus[1])
  rownames(p_n_l) <- paste0("n_r=", nvec)
  colnames(p_n_l) <- paste0("n_l=", nvec)
  p_r_cs_l_cp <-
    p_n_r * p_n_l * (n_r / (n_r + n_l)) * (1 - (one_minus_p_b_1) ^ n_l)
  p_l <- sum(p_r_cs_l_cp) # branch 2 survives till the present
  p_r_cp_l_cp <- p_r_cs_l_cp * (1 - (one_minus_p_b_1) ^ (n_r - 1))
  p_rl  <- sum(p_r_cp_l_cp) # both branches 1 and 2 survive till the present

  pc_1 <- 2 * p_s * p_l + 2 * (1 - p_s) * p_rl
  pc_2 <- 2 * p_s * p_l
  pc_3 <- 2 * p_s * p_rl

  pc <- (cond == 0) * 1 +
    (cond == 1) * pc_1 +
    (cond == 2) * pc_2 +
    (cond == 3) * pc_3

  if (pc == 0) {
    print("Conditional probability is equal to zero!")
  }
  # there might be small numerical errors
  if (pc > 1 && pc < (1 + 1e13)) {
    pc <- 1
  }
  testit::assert(pc > 0, fact = paste0("Pc <= 0 for parameters ", paste(c(pars_m, pars_s), collapse = " ")))
  testit::assert(pc <= 1, fact = paste0("Pc > 1 for parameters ", paste(c(pars_m, pars_s), collapse = " ")))
  return(pc)
}
