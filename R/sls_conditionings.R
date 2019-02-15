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
  t_d <- brts_s1[1]
  t_c <- brts_m1[1]
  t_p <- 0
  aa <- abs(t_d - t_c); bb <- abs(t_p - t_d)

  if (n_0 != 2) {
    stop("Pc can be calculated only if phylogeny starts with a crown!")
  }

  p_s <- sls::p_t(
    t = bb,
    lambda = lambdas[2],
    mu = mus[2]
  )

  nvec <- 1:n_max
  ns1  <- row(matrix(NA, nrow = n_max, ncol = n_max))
  ns2  <- col(matrix(NA, nrow = n_max, ncol = n_max))
  p_a   <- sls::p_t(t = aa, lambda = lambdas[1], mu = mus[1]); p_a
  u_a   <- sls::ut(t = aa, lambda = lambdas[1], mu = mus[1]); u_a
  p_b1  <- sls::p_t(t = bb, lambda = lambdas[1], mu = mus[1]); p_b1
  one_minus_p_b1  <- sls::one_minus_pt(t = bb, lambda = lambdas[1], mu = mus[1])
  p_b2  <- sls::p_t(t = bb, lambda = lambdas[2], mu = mus[2]); p_b2
  p_ns1 <- sls::pn(n = ns1, t = aa, lambda = lambdas[1], mu = mus[1])
  rownames(p_ns1) <- paste0("ns1=", nvec)
  colnames(p_ns1) <- paste0("ns2=", nvec)
  p_ns2 <- sls::pn(n = ns2, t = aa, lambda = lambdas[1], mu = mus[1])
  rownames(p_ns2) <- paste0("ns1=", nvec)
  colnames(p_ns2) <- paste0("ns2=", nvec)
  aux1 <- p_ns1 * p_ns2 * (ns1 / (ns1 + ns2)) * (1 - (one_minus_p_b1) ^ ns2)
  p_1   <- sum(aux1) #branch 2 survives till the present
  aux2 <- aux1 * (1 - (one_minus_p_b1) ^ (ns1 - 1))
  p_2   <- sum(aux2) #both branches 1 and 2 survive till the present

  pc_1  <- 2 * p_s * p_1 + 2 * (1 - p_s) * p_2
  pc_2  <- 2 * p_s * p_1
  pc_3  <- 2 * p_s * p_2

  pc <- (cond == 0) * 1 +
    (cond == 1) * pc_1 +
    (cond == 2) * pc_2 +
    (cond == 3) * pc_3

  if (pc == 0) {
    print("Conditional probability is equal to zero!")
  }
  testit::assert(pc > 0)
  testit::assert(pc <= 1)

  return(pc)
}
