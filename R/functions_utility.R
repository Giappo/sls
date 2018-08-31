#' @title Brts to time intervals
#' @author Giovanni Laudanno
#' @description Calculates time intervals given the branching times
#' @inheritParams default_params_doc
#' @return The time intervals
#' @export
brts2time_intervals = function (brts){
  time_points <- -unlist(unname(sort(abs(brts), decreasing = TRUE)) )
  branching_times <- -sort(abs(as.numeric(time_points)),decreasing = TRUE)
  births <- c(0, unname(table(branching_times))[-1])
  unique_branching_times <- as.numeric(names(table(branching_times)))
  time_intervals <- c((diff(unique_branching_times)), (abs(tail(unique_branching_times,1))) )
  births <- births[-1]
  time_intervals <- c(0, time_intervals)
  births <- c(0, births)
  return(time_intervals)
}

#' @title Transition matrix builder
#' @author Giovanni Laudanno
#' @description Builds the transition matrix to integrate the differential equations of the P-equation
#' @inheritParams default_params_doc
#' @return The transition matrix
#' @export
P_transition_matrix <- function(lambda, mu, matrix_size){
  nvec <- 0:(matrix_size - 1)
  M <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  M[row(M) == col(M) + 1] <- M[row(M) == col(M) + 1] + lambda * (nvec[1:(matrix_size - 1)])
  M[row(M) == col(M) - 1] <- M[row(M) == col(M) - 1] + mu * (nvec[2:matrix_size])
  M[row(M) == col(M)] <- M[row(M) == col(M)] - (lambda + mu) * (nvec[1:matrix_size])
  M[matrix_size, matrix_size] <- - M[matrix_size - 1, matrix_size]; M
  testit::assert(colSums(M) < 1e-10)
  return(M)
}

#' @title L to Lineage-Through-Time
#' @author Giovanni Laudanno
#' @description From an L matrix it yields the corresponding lineage-through-time function.
#' @inheritParams default_params_doc
#' @return The function returns a matrix:
#' \itemize{
#' \item the first row contains the time points when the amount of lineages change;
#' \item the second row gives the number of lineages immediately after that time point;
#' }
#' @export
L2LTT <- function(L, shifts = NULL) {
  #full LTT
  adds <- L[, 1]
  rems <- L[L[,4] != -1, 4]
  times <- sort(c(adds, rems), decreasing = TRUE)
  kvec <- rep(NA, length(times))
  if (length(kvec) > 0)
  {
    kvec[1] <- 1
    for (j in seq_along(times)[-1])
    {
      kvec[j] <- kvec[j - 1] + times[j] %in% adds - times[j] %in% rems
    }; kvec
    kmat <- rbind(times, kvec)
    if (sum(kmat[1,1] == kmat[1,]) > 1) {kmat <- kmat[, -(1:(sum(kmat[1,1] == kmat[1,]) - 1))]}; kmat
  }else
  {
    kmat <- kvec
  }
  LTT_full <- kmat
  dim(LTT_full) <- c(2, length(unique(times)))

  #reconstructed LTT
  extant   <- which(L[, 4] == -1)
  Lextant  <- L[extant ,]; dim(Lextant)  <- c(length(extant) , ncol(L))
  if (!is.null(shifts))
  {
    shifted <- which(L[, 4] %in% abs(shifts))
    Lshifted <- L[shifted,]; dim(Lshifted) <- c(length(shifted), ncol(L))
    L1 <- rbind(Lextant, Lshifted)
  }else
  {
    L1 <- Lextant
  }
  L1 <- L1[order(abs(L1[, 3])),]
  adds1 <- L1[, 1]
  rems1 <- L1[L1[,4] != -1, 4]
  times1 <- sort(c(adds1, rems1), decreasing = TRUE)
  kvec1 <- rep(NA, length(times1))
  if (length(kvec1) > 0)
  {
    kvec1[1] <- 1
    for (j in seq_along(times1)[-1])
    {
      kvec1[j] <- kvec1[j - 1] + times1[j] %in% adds1 - times1[j] %in% rems1
    }; kvec1
    kmat1 <- rbind(times1, kvec1)
    if (sum(kmat1[1,1] == kmat1[1,]) > 1) {kmat1 <- kmat1[, -(1:(sum(kmat1[1,1] == kmat1[1,]) - 1))]}; kmat1
  }else
  {
   kmat1 <- kvec1
  }
  LTT_reconstructed <- kmat1; reco_LTT; L1
  dim(LTT_reconstructed) <- c(2, length(unique(times1)))

  return(list(LTT_full = LTT_full, LTT_reconstructed = LTT_reconstructed))
}

#####OTHER FUNCTIONS THAT I DON'T NEED/USE NOW
#' #' Does something
#' #' @inheritParams default_params_doc
#' #' @return result
#' #' @export
#' old_Pc_1shift <- function(brtsM, brtsS, lambdas, mus) {
#'
#'   tp <- 0 ;tc <- brtsM[1]; ts <- brtsS[1]
#'   testit:::assert(tp > ts)
#'   testit:::assert(ts > tc)
#'
#'   PM1   <- (1 - sls::pn(n = 0, t = tp - tc, lambda = lambdas[1], mu = mus[1]))
#'   PM2cs <- (1 - sls::pn(n = 0, t = ts - tc, lambda = lambdas[1], mu = mus[1]))
#'   PS    <- (1 - sls::pn(n = 0, t = tp - ts, lambda = lambdas[2], mu = mus[2]))
#'
#'   P1cs <- sls::pn(n = 1, t = ts - tc, lambda = lambdas[1], mu = mus[1])
#'   P0sp <- sls::pn(n = 0, t = tp - ts, lambda = lambdas[1], mu = mus[1])
#'   ucs  <- sls::ut(t = ts - tc, lambda = lambdas[1], mu = mus[1])
#'
#'   PM2cp <- P1cs *
#'     (1 - P0sp) *
#'     ucs *
#'     (1 - ucs)^-1 *
#'     (1 - ucs * P0sp)^-1
#'
#'   Pc1 <- PM1 * (PS * PM2cs + (1 - PS) * PM2cp)
#'
#'   Pc2 <- PM1 * PS * PM2cp
#'   Pc3 <- PM1 * PS * PM2cs
#'   return(list(Pc1 = Pc1, Pc2 = Pc2, Pc3 = Pc3))
#' }
#'
#' #' Does something
#' #' @inheritParams default_params_doc
#' #' @return result
#' #' @export
#' Bart_Pc1 <- function(brtsM, brtsS, lambdas, mus) {
#'   tc <- brtsM[1]
#'   ts <- brtsS[1]
#'   tp <- 0
#'   A <- ts - tc; B <- tp - ts
#'   Nmax <- 100
#'   nvec <- 1:Nmax
#'   ns1  <- row(matrix(NA, nrow = Nmax, ncol = Nmax))
#'   ns2  <- col(matrix(NA, nrow = Nmax, ncol = Nmax))
#'   pA   <- sls::pt(t = A, lambda = lambdas[1], mu = mus[1]); pA
#'   uA   <- sls::ut(t = A, lambda = lambdas[1], mu = mus[1]); uA
#'   pB   <- sls::pt(t = B, lambda = lambdas[1], mu = mus[1]); pB
#'   pns1 <- sls::pn(n = ns1, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns1) <- paste0("ns1=", nvec); colnames(pns1) <- paste0("ns2=", nvec); pns1
#'   pns2 <- sls::pn(n = ns2, t = A, lambda = lambdas[1], mu = mus[1]); rownames(pns2) <- paste0("ns1=", nvec); colnames(pns2) <- paste0("ns2=", nvec); pns2
#'   aux1 <- pns1 * pns2 * (ns1/(ns1 + ns2)) * (1 - (1 - pB)^ns2)
#'   aux2 <- pns1 * pns2 * (ns2/(ns1 + ns2)) * (1 - (1 - pB)^ns1)
#'   PS <- 1 - pn(n = 0, t = B, lambda = lambdas[2], mu = mus[2])
#'   out  <- PS * (sum(aux1) + sum(aux2)); out
#'   return(out)
#' }
#' #' Does something
#' #' @inheritParams default_params_doc
#' #' @return result
#' #' @export
#' Pc_1shift2 <- function(brtsM, brtsS, lambdas, mus) {
#'
#'   tp <- 0 ;tc <- brtsM[1]; ts <- brtsS[1]
#'   testit:::assert(tp > ts)
#'   testit:::assert(ts > tc)
#'
#'   PM1   <- (1 - sls::pn(n = 0, t = tp - tc, lambda = lambdas[1], mu = mus[1]))
#'   PM2cs <- (1 - sls::pn(n = 0, t = ts - tc, lambda = lambdas[1], mu = mus[1]))
#'   PS    <- (1 - sls::pn(n = 0, t = tp - ts, lambda = lambdas[2], mu = mus[2]))
#'
#'   P1cs <- sls::pn(n = 1, t = ts - tc, lambda = lambdas[1], mu = mus[1])
#'   P0sp <- sls::pn(n = 0, t = tp - ts, lambda = lambdas[1], mu = mus[1])
#'   ucs  <- sls::ut(t = ts - tc, lambda = lambdas[1], mu = mus[1])
#'
#'   PM2cp <- P1cs *
#'     (1 - P0sp) *
#'     ucs *
#'     (1 - ucs)^-1 *
#'     (1 - ucs * P0sp)^-1; PM2cp
#'
#'   PM2cp2 <- P1cs *(
#'     (1 - ucs)^-2 -
#'     (1 - ucs * P0sp)^-2); PM2cp2
#'
#'   # # ######### summing - divide by ns
#'   # ns <- 2:1e5; PM2cp3 <- P1cs *
#'   #   sum(
#'   #     ucs^(ns - 1) *
#'   #     (1 - P0sp^(ns - 1))#do i divide or not? In our framework I divide and simplify this
#'   #   ); PM2cp3
#'   # # ###########
#'   #
#'   # # ######### summing - don't divide by ns
#'   # ns <- 2:1e5; PM2cp4 <- P1cs *
#'   #   sum(
#'   #     ucs^(ns - 1) *
#'   #       (1 - P0sp^(ns - 1)) * #do i divide or not? In our framework I divide and simplify this
#'   #     ns
#'   #   ); PM2cp4
#'   # # ###########
#'
#'   Pc1 <- PM1 * (PS * PM2cs + (1 - PS) * PM2cp2)
#'   Pc2 <- PM1 * PS * PM2cp2
#'   Pc3 <- PM1 * PS * PM2cs
#'   return(list(Pc1 = Pc1, Pc2 = Pc2, Pc3 = Pc3))
#' }
#'
#' #' Does something
#' #' @inheritParams default_params_doc
#' #' @return result
#' #' @export
#' Qc3_1shift0 <- function(brtsM, brtsS, lambdas, mus, lx = 500) {
#'
#'   #BASIC SETTINGS AND CHECKS
#'   Nclades <- length(lambdas)
#'   brts_list <- list(brtsM = rep(brtsM[1], sum(brtsM == brtsM[1])), brtsS = brtsS[1])
#'   shift_times <- unlist(lapply(brts_list, FUN = function(x) x[1]))
#'   abstol <- 1e-16; reltol <- 1e-10
#'   #MAIN
#'
#'   #ADJUSTING DATA
#'   nvec <- 0:lx
#'   clade <- 0 #clade == 1 is the main clade, clade == 2 is the subclade
#'   logliks <- rep(NA, Nclades)
#'   Pm <- CC <- DD <- list(0, Nclades)
#'   #LIKELIHOOD INTEGRATION
#'   while ((clade <- clade + 1) <= Nclades)
#'   {
#'     #SETTING CLADE CONDITIONS
#'     shift_times2 <- shift_times[shift_times > min(brts_list[[clade]])]
#'     time_points <- sort(unique(c(brts_list[[clade]], shift_times2)), decreasing = FALSE)
#'     time_intervals <- sls::brts2time_intervals(time_points)
#'     lambda <- lambdas[clade]
#'     mu     <- mus[clade]
#'     N0     <- sum(brts_list[[clade]] == brts_list[[clade]][1])
#'
#'     #SETTING INITIAL CONDITIONS (there's always a +1 because of Q0)
#'     Qi <- c(1, rep(0, lx))
#'     Qt <- matrix(0, ncol = (lx + 1), nrow = length(time_intervals))
#'     Qt[1,] <- Qi
#'     dimnames(Qt)[[2]] <- paste0("Q", 0:lx)
#'     k <- N0
#'     t <- 2
#'     D <- C <- rep(1, (length(time_intervals)))
#'
#'     #EVOLVING THE INITIAL STATE TO THE LAST BRANCHING POINT
#'     while (t <= length(time_intervals))
#'     {
#'       #Applying A operator
#'       if (lambda == 0 && mu == 0)
#'       {
#'         Qt[t,] <- Qt[(t-1),]
#'       }else
#'       {
#'         transition_matrix <- DDD:::dd_loglik_M_aux(pars = c(lambda, mu, Inf), lx = lx + 1, k = k, ddep = 1)
#'         Qt[t,] <- abs(expoRkit::expv(v = Qt[(t-1),], x = transition_matrix, t = time_intervals[t]))
#'       }
#'       # if (methode == "analytical") {
#'       # }
#'       # else {
#'       #   transition_matrix <- MBD:::create_A(lambda = lambda, mu = mu, nu = 0, q = 0, k = k, max_number_of_species = lx)
#'       #   Qt[t,] <- deSolve::ode(y = Qt[(t-1),], times = c(0, time_intervals[t]), func = MBD:::mbd_loglik_rhs,
#'       #                          parms = transition_matrix, atol = 1e-16, rtol = 1e-10)[2,-1]
#'       # }
#'
#'       #Applying C operator (this is a trick to avoid precision issues)
#'       C[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * C[t]
#'
#'       #what time is it?
#'       tempo <- time_points[t]
#'
#'       if (t < length(time_intervals))
#'       {
#'         #Applying B operator
#'         if (all(tempo != shift_times))
#'         {
#'           Qt[t,] <- Qt[t,] * k * lambda
#'           k <- k + 1
#'         }else
#'         {
#'           # Qt[t,] <- Qt[t,] * (k + nvec)^-1 #first version
#'           # Qt[t,] <- Qt[t,] * nvec * (k + nvec)^-1 #second version
#'           Qt[t,] <- Qt[t,] * k * (k + nvec)^-1 #third version
#'           k <- k - 1
#'         }
#'
#'         #Applying D operator (this works exactly like C)
#'         D[t] <- 1/(sum(Qt[t,])); Qt[t,] <- Qt[t,] * D[t]
#'
#'         #Updating running parameters
#'         t <- t + 1
#'       }else{break}
#'     }
#'
#'     #Selecting the state I am interested in
#'     vm <- 1/choose((k + nvec), k)
#'     Pm[[clade]] <- vm * Qt[t, (nvec + 1)] #I have to include +1 because of Q0
#'
#'     #Removing C and D effects from the LL
#'     # loglik <- log(P) - sum(log(C)) - sum(log(D))
#'     CC[[clade]] <- C
#'     DD[[clade]] <- D
#'
#'     #Various checks
#'     # loglik <- as.numeric(loglik)
#'     # if (is.nan(loglik) | is.na(loglik))
#'     # {
#'     #   loglik <- -Inf
#'     # }
#'     # logliks[clade] <- loglik
#'   }
#'   # total_loglik <- sum(logliks)
#'
#'   Pc <- sum(Pm[[1]]) * sum(CC[[1]])^-1 * sum(DD[[1]])^-1 +
#'     sum(Pm[[2]]) * sum(CC[[2]])^-1 * sum(DD[[2]])^-1
#'
#'   # Pc <- 2 * (PM.1.12 + PS * PM.1.2)
#'   # return(total_loglik)
#'   return(Pc)
#' }
#
#' #' Does something
#' #' @inheritParams default_params_doc
#' #' @return result
#' #' @export
#' Qc3_1shift <- function(brtsM, brtsS, lambdas, mus, maxN = 10) {
#'
#'   brts_list <- list(brtsM = rep(brtsM[1], sum(brtsM == brtsM[1])), brtsS = brtsS[1])
#'   shift_times <- unlist(lapply(brts_list, FUN = function(x) x[1]))
#'   abstol <- 1e-16; reltol <- 1e-10
#'   tc <- brts_list$brtsM[1]; ts <- brts_list$brtsS[1]; tp <- 0
#'   lambda <- lambdas[1]; mu <- mus[1]
#'
#'   nvec1 <- 0:maxN; maxN1 <- length(nvec1); maxN1
#'   nvec2 <- 0:((maxN2 <- maxN1^2) - 1); maxN2
#'   N0 <- sum(brts_list[[1]] == brts_list[[1]][1]); N0
#'
#'   Pi <- c(rep(0, maxN2)); Pi[N0 + 1] <- 1
#'   Pt <- matrix(0, ncol = maxN2, nrow = 2)
#'   Pt[1,] <- Pi
#'   dimnames(Pt)[[2]] <- paste0("P", nvec2); head(Pt)
#'   t <- 2
#'   D <- C <- rep(1, 2)
#'
#'   transition_matrix <- sls::P_transition_matrix(lambda = lambda, mu = mu, matrix_size = maxN2); dim(transition_matrix)
#'   testit::assert(dim(transition_matrix) == ncol(Pt))
#'   testit::assert(is.na(transition_matrix) == 0)
#'   times <- c(0, ts - tc)
#'   Pt[t,] <- abs(expoRkit::expv(v = Pt[(t-1),], x = transition_matrix, t = (ts - tc)))
#'   grid_n <- expand.grid(n1 = nvec1, n2 = nvec1); dim(grid_n)
#'   # grid_n <- grid_n[apply(grid_n, MARGIN = 1, FUN = function(x) sum(x) < max(nvec)),]
#'   Pplus <- Pt[t, grid_n[,1] + grid_n[,2] + 1 + 1] *
#'            (grid_n[,1] + 1)/((grid_n[,1] + grid_n[,2]) * (grid_n[,1] + grid_n[,2] + 1))
#'   Pplus[1] <- 0
#'   names(Pplus) <- paste0("P(", grid_n[,1],",", grid_n[,2],")")
#'
#'   times <- c(0, tp - ts)
#'   Pf <- deSolve::ode(y = Pplus, times = times, func = sls::Pn1n2_rhs,
#'   parms = c(lambda, mu), atol = 1e-16, rtol = 1e-10)[2,-1]; length(Pf)
#'
#'   Pfmatrix <- matrix(Pf, nrow = maxN1, ncol = maxN1)
#'   colnames(Pfmatrix) <- paste0("n2=", nvec1); rownames(Pfmatrix) <- paste0("n1=", nvec1)
#'   PM.1.2  <- sum(Pfmatrix[1,-1])
#'   PM.1.12 <- sum(Pfmatrix[-1,-1])
#'   PS <- (1 - sls::pn(n = 0, t = tp - ts, lambda = lambdas[2], mu = mus[2]))
#'
#'   Pc <- 2 * (PM.1.12 + PS * PM.1.2)
#'
#'   return(Pc)
#' }
#'
#' #' Does something
#' #' @inheritParams default_params_doc
#' #' @return result
#' #' @export
#' Pn1n2_rhs <- function(t, x, pars) {
#'
#'   with(as.list(x), {
#'     maxN <- sqrt(length(x)) - 1; lambda <- pars[1]; mu <- pars[2]
#'     Gn <- expand.grid(n1 = 0:maxN, n2 = 0:maxN); dim(Gn) # Gn <- grid_n[apply(grid_n, MARGIN = 1, FUN = function(x) sum(x) < max(nvec)),]
#'
#'     i0 <- 1:nrow(Gn)
#'     i_n1_minus <- i0 + rep(c(0, rep(-1, maxN)), maxN + 1)
#'     i_n1_plus  <- i0 + rep(c(rep(+1, maxN), 0), maxN + 1)
#'     i_n2_minus <- i0 + c(rep(0, maxN + 1), rep(rep(-(maxN + 1), maxN + 1), maxN))
#'     i_n2_plus  <- i0 + c(rep(rep(+(maxN + 1), maxN + 1), maxN), rep(0, maxN + 1))
#'
#'     n1      <- Gn[i0, 1]
#'     n2      <- Gn[i0, 2]
#'     n1minus <- Gn[i_n1_minus, 1]
#'     n1plus  <- Gn[i_n1_plus, 1]
#'     n2minus <- Gn[i_n2_minus, 2]
#'     n2plus  <- Gn[i_n2_plus, 2]
#'
#'     testit::assert(n1minus[n1 != 0] - n1[n1 != 0] == -1)
#'     testit::assert(n1plus[n1 != maxN] - n1[n1 != maxN] == 1)
#'     testit::assert(n2minus[n2 != 0] - n2[n2 != 0] == -1)
#'     testit::assert(n2plus[n2 != maxN] - n2[n2 != maxN] == 1)
#'
#'     testit::assert(length(i_n1_plus)  == length(i0))
#'     testit::assert(length(i_n2_plus)  == length(i0))
#'     testit::assert(length(i_n1_minus) == length(i0))
#'     testit::assert(length(i_n2_minus) == length(i0))
#'     testit::assert(length(n1     ) == length(i0))
#'     testit::assert(length(n2     ) == length(i0))
#'     testit::assert(length(n1minus) == length(i0))
#'     testit::assert(length(n2minus) == length(i0))
#'     testit::assert(length(n1plus ) == length(i0))
#'     testit::assert(length(n2plus ) == length(i0))
#'     testit::assert(Gn[i_n1_minus[Gn[i0, 1] != 0], 1] - Gn[i0[Gn[i0, 1] != 0], 1] == -1)
#'     testit::assert(Gn[i_n1_plus[Gn[i0, 1] != maxN], 1] - Gn[i0[Gn[i0, 1] != maxN], 1] == 1)
#'     testit::assert(Gn[i_n2_minus[Gn[i0, 2] != 0], 2] - Gn[i0[Gn[i0, 2] != 0], 2] == -1)
#'     testit::assert(Gn[i_n2_plus[Gn[i0, 2] != maxN], 2] - Gn[i0[Gn[i0, 2] != maxN], 2] == 1)
#'     lavec <- rep(1, 2 * maxN^2) * lambda
#'     muvec <- rep(1, 2 * maxN^2) * mu
#'
#'     dx = lavec[i0] * (n1minus * x[i_n1_minus] + n2minus * x[i_n2_minus]) +
#'       muvec[i0] * (n1plus * x[i_n1_plus] + n2plus * x[i_n2_plus]) +
#'       -(muvec[i0] + lavec[i0]) * (n1 + n2) * x[i0]
#'
#'     testit::assert(sum(dx) < (lambda + mu) * 1e-2)
#'     names(dx) = names(x)
#'     return(list(dx))
#'   })
#' }
#
# GG <- function(z, t, lambda, mu) {
#   LL  <- exp((mu - lambda) * t)
#   out <- (mu - mu * z - LL * (mu - lambda * z)) *
#          (lambda - lambda * z - LL * (mu - lambda * z))^-1
#   return(out)
# }
# Pc1_1shift <- function(brtsM, brtsS, lambdas, mus) {
#
#   tp <- 0 ;tc <- brtsM[1]; ts <- brtsS[1]
#   testit:::assert(tp > ts)
#   testit:::assert(ts > tc)
#
#   PM1   <- (1 - sls::pn(n = 0, t = tp - tc, lambda = lambdas[1], mu = mus[1]))
#   PM2cs <- (1 - sls::pn(n = 0, t = ts - tc, lambda = lambdas[1], mu = mus[1]))
#   PS    <- (1 - sls::pn(n = 0, t = tp - ts, lambda = lambdas[2], mu = mus[2]))
#
#   P1cs <- sls::pn(n = 1, t = ts - tc, lambda = lambdas[1], mu = mus[1])
#   P0sp <- sls::pn(n = 0, t = tp - ts, lambda = lambdas[1], mu = mus[1])
#   ucs  <- sls::ut(t = ts - tc, lambda = lambdas[1], mu = mus[1])
#
#   PM2cp <- P1cs *
#     (1 - P0sp) *
#     ucs *
#     (1 - ucs)^-1 *
#     (1 - ucs * P0sp)^-1
#
#   out <- PM1 * (PS * PM2cs + (1 - PS) * PM2cp)
#   return(out)
# }
# Pc3_1shift <- function(brtsM, brtsS, lambdas, mus) {
#
#   tp <- 0 ;tc <- brtsM[1]; ts <- brtsS[1]
#   testit:::assert(tp > ts)
#   testit:::assert(ts > tc)
#   Pc <- (1 - sls::pn(n = 0, t = tp - tc, lambda = lambdas[1], mu = mus[1])) * #survival of crown species not undergoing the shift
#     (1 - sls::pn(n = 0, t = ts - tc, lambda = lambdas[1], mu = mus[1])) * #survival of mainclade til shift time
#     (1 - sls::pn(n = 0, t = tp - ts, lambda = lambdas[2], mu = mus[2]))   #survival of subclade from shift time to present
#   return(Pc)
# }
# Pn1n2_matrix <- function(lambda, mu, matrix_size) {
#   #so far it's just a copy of Pn1n2_rhs. Still needs some manipulations
#   maxN <- sqrt(matrix_size) - 1; lambda <- pars[1]; mu <- pars[2]
#   Gn <- expand.grid(n1 = 0:maxN, n2 = 0:maxN); dim(Gn) # Gn <- grid_n[apply(grid_n, MARGIN = 1, FUN = function(x) sum(x) < max(nvec)),]
#
#   i0 <- 1:nrow(Gn)
#   i_n1_minus <- i0 + rep(c(0, rep(-1, maxN)), maxN + 1)
#   i_n1_plus  <- i0 + rep(c(rep(+1, maxN), 0), maxN + 1)
#   i_n2_minus <- i0 + c(rep(0, maxN + 1), rep(rep(-(maxN + 1), maxN + 1), maxN))
#   i_n2_plus  <- i0 + c(rep(rep(+(maxN + 1), maxN + 1), maxN), rep(0, maxN + 1))
#
#   n1      <- Gn[i0, 1]
#   n2      <- Gn[i0, 2]
#   n1minus <- Gn[i_n1_minus, 1]
#   n1plus  <- Gn[i_n1_plus, 1]
#   n2minus <- Gn[i_n2_minus, 2]
#   n2plus  <- Gn[i_n2_plus, 2]
#
#   testit::assert(n1minus[n1 != 0] - n1[n1 != 0] == -1)
#   testit::assert(n1plus[n1 != maxN] - n1[n1 != maxN] == 1)
#   testit::assert(n2minus[n2 != 0] - n2[n2 != 0] == -1)
#   testit::assert(n2plus[n2 != maxN] - n2[n2 != maxN] == 1)
#
#   testit::assert(length(i_n1_plus)  == length(i0))
#   testit::assert(length(i_n2_plus)  == length(i0))
#   testit::assert(length(i_n1_minus) == length(i0))
#   testit::assert(length(i_n2_minus) == length(i0))
#   testit::assert(length(n1     ) == length(i0))
#   testit::assert(length(n2     ) == length(i0))
#   testit::assert(length(n1minus) == length(i0))
#   testit::assert(length(n2minus) == length(i0))
#   testit::assert(length(n1plus ) == length(i0))
#   testit::assert(length(n2plus ) == length(i0))
#   testit::assert(Gn[i_n1_minus[Gn[i0, 1] != 0], 1] - Gn[i0[Gn[i0, 1] != 0], 1] == -1)
#   testit::assert(Gn[i_n1_plus[Gn[i0, 1] != maxN], 1] - Gn[i0[Gn[i0, 1] != maxN], 1] == 1)
#   testit::assert(Gn[i_n2_minus[Gn[i0, 2] != 0], 2] - Gn[i0[Gn[i0, 2] != 0], 2] == -1)
#   testit::assert(Gn[i_n2_plus[Gn[i0, 2] != maxN], 2] - Gn[i0[Gn[i0, 2] != maxN], 2] == 1)
#   lavec <- rep(1, 2 * maxN^2) * lambda
#   muvec <- rep(1, 2 * maxN^2) * mu
#
#   Gn[2,]
#
#   dx = lavec[i0] * (n1minus * x[i_n1_minus] + n2minus * x[i_n2_minus]) +
#     muvec[i0] * (n1plus  * x[i_n1_plus ] + n2plus  * x[i_n2_plus ]) +
#     -(muvec[i0] + lavec[i0]) * (n1 + n2) * x[i0]
#
#   testit::assert(sum(dx) < 1e-3)
#   names(dx) = names(x)
# }

