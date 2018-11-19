#' @title sls-conditioning
#' @author Giovanni Laudanno
#' @description Calculates three different kind of conditioning probabilities
#' @inheritParams default_params_doc
#' @return Three different conditioning probabilities
#' @export
pc_1shift_old <- function(pars1, pars2, brtsM, brtsS) {

  missnumspec <- c(0,0)
  lambdas <- c(pars1[1], pars1[4])
  mus     <- c(pars1[2], pars1[5])
  Ks      <- c(pars1[3], pars1[6])
  td      <- pars1[7]
  td      <- abs(td) * sign(brtsM[1])

  lx     <- pars2[1]
  cond   <- pars2[3]
  tsplit <- pars2[4]
  soc    <- pars2[6]

  tc <- brtsM[1]
  ts <- td
  tp <- 0
  A <- abs(ts - tc); B <- abs(tp - ts)

  if (soc != 2) {stop("Pc can be calculated only if phylogeny starts with a crown!")}

  PS   <- 1 - pn(n = 0, t = B, lambda = lambdas[2], mu = mus[2])
  nvec <- 1:lx
  ns1  <- row(matrix(NA, nrow = lx, ncol = lx))
  ns2  <- col(matrix(NA, nrow = lx, ncol = lx))
  pA   <- sls::pt(t = A, lambda = lambdas[1], mu = mus[1]); pA
  uA   <- sls::ut(t = A, lambda = lambdas[1], mu = mus[1]); uA
  pB1  <- sls::pt(t = B, lambda = lambdas[1], mu = mus[1]); pB1
  pB2  <- sls::pt(t = B, lambda = lambdas[2], mu = mus[2]); pB2
  pns1 <- sls::pn(n = ns1, t = A, lambda = lambdas[1], mu = mus[1])
  rownames(pns1) <- paste0("ns1=", nvec);
  colnames(pns1) <- paste0("ns2=", nvec); head(pns1)
  pns2 <- sls::pn(n = ns2, t = A, lambda = lambdas[1], mu = mus[1])
  rownames(pns2) <- paste0("ns1=", nvec);
  colnames(pns2) <- paste0("ns2=", nvec); head(pns2)
  aux1 <- pns1 * pns2 * (ns1 / (ns1 + ns2)) * (1 - (1 - pB1)^ns2)
  P1   <- sum(aux1) #branch 2 survives till the present
  aux2 <- aux1 * (1 - (1 - pB1) ^ (ns1 - 1)); head(aux2)
  P2   <- sum(aux2) #both branches 1 and 2 survive till the present

  pc_1  <- 2 * PS * P1 + 2 * (1 - PS) * P2
  pc_2  <- 2 * PS * P2
  pc_3  <- 2 * PS * P1

  return(list(pc_1 = pc_1, pc_2 = pc_2, pc_3 = pc_3))
}
