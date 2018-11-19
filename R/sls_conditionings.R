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
  nmax = 1e2,
  n_0 = 2
) {
  lambdas <- c(pars_m[1], pars_s[1])
  mus     <- c(pars_m[2], pars_s[2])

  brts_m1 <- sort(abs(brts_m), decreasing = TRUE)
  brts_s1 <- sort(abs(brts_s), decreasing = TRUE)
  td <- brts_s1[1]

  tc <- brts_m1[1]
  tp <- 0
  A <- abs(td - tc); B <- abs(tp - td)

  if (n_0 != 2) {
    stop("Pc can be calculated only if phylogeny starts with a crown!")
  }

  PS   <- 1 - sls::pn(
    n = 0,
    t = B,
    lambda = lambdas[2],
    mu = mus[2]
  )
  nvec <- 1:nmax
  ns1  <- row(matrix(NA, nrow = nmax, ncol = nmax))
  ns2  <- col(matrix(NA, nrow = nmax, ncol = nmax))
  p_a   <- sls::pt(t = A, lambda = lambdas[1], mu = mus[1]); p_a
  u_a   <- sls::ut(t = A, lambda = lambdas[1], mu = mus[1]); u_a
  p_b1  <- sls::pt(t = B, lambda = lambdas[1], mu = mus[1]); p_b1
  p_b2  <- sls::pt(t = B, lambda = lambdas[2], mu = mus[2]); p_b2
  p_ns1 <- sls::pn(n = ns1, t = A, lambda = lambdas[1], mu = mus[1])
  rownames(p_ns1) <- paste0("ns1=", nvec)
  colnames(p_ns1) <- paste0("ns2=", nvec)
  p_ns2 <- sls::pn(n = ns2, t = A, lambda = lambdas[1], mu = mus[1])
  rownames(p_ns2) <- paste0("ns1=", nvec)
  colnames(p_ns2) <- paste0("ns2=", nvec)
  aux1 <- p_ns1 * p_ns2 * (ns1 / (ns1 + ns2)) * (1 - (1 - p_b1) ^ ns2)
  P1   <- sum(aux1) #branch 2 survives till the present
  aux2 <- aux1 * (1 - (1 - p_b1) ^ (ns1 - 1)); head(aux2)
  P2   <- sum(aux2) #both branches 1 and 2 survive till the present

  pc_1  <- 2 * PS * P1 + 2 * (1 - PS) * P2
  pc_4  <- 2 * PS * P2
  pc_3  <- 2 * PS * P1

  pc <- (cond == 0) * 1 +
    (cond == 1) * pc_1 +
    (cond == 3) * pc_3 +
    (cond == 4) * pc_4

  return(pc)
}
