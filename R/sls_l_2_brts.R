#' @title Convert l_0 to brts
#' @author Giovanni Laudanno
#' @description Convert l_0 to brts
#' @inheritParams default_params_doc
#' @return branching times
#' @export
sls_l_2_brts <- function(
  l_0,
  n_0,
  dropextinct = TRUE
) {

  n_row_l <- nrow(l_0)
  n_col_l <- ncol(l_0)
  brts <- NULL
  l_0 <- l_0[order(abs(l_0[, 3])), 1:n_col_l]
  dim(l_0) <- c(n_row_l, n_col_l)
  age <- l_0[1, 1]
  l_0[, 1] <- age - l_0[, 1]
  l_0[1, 1] <- -1
  not_min_1 <- l_0[, 4] != -1 & l_0[, 5] == 0
  extant <- abs(l_0[, 3])[!not_min_1]
  if (any(l_0[, 5] > 0)) {
    shifted <- which(l_0[, 5] > 0)
    shift_times <- rep(0, length(extant))
    shift_times[extant == shifted] <- l_0[shifted, 4]
  } else {
    shifted <- 0
    shift_times <- rep(0, length(extant))
  }

  l_0[not_min_1, 4] <- age - l_0[not_min_1, 4]
  if (dropextinct == TRUE) {
    sall <- abs(l_0[extant, 3])
    t_end <- age * (sall != shifted) + shift_times * (sall == shifted)
  } else {
    sall <- which(l_0[, 4] >= -1)
    t_end <- (l_0[, 4] == -1) * age + (l_0[, 4] > -1) * l_0[, 4]
  }
  l_0_redux <- l_0[, -4]
  dim(l_0_redux) <- c(n_row_l, n_col_l - 1)
  df <- data.frame(
    matrix(
      l_0_redux[sall, ],
      nrow = length(sall),
      ncol = n_col_l - 1
    )
  )
  lin_list <- cbind(
    df,
    paste0("t", abs(l_0_redux[sall, 3])),
    t_end
  )
  names(lin_list) <- c(
    "birth",
    "parent",
    "id",
    "shift",
    "label",
    "t_end"
  )
  lin_list$label <- as.character(lin_list$label)
  if (nrow(lin_list) == 1) {
    brts <- age
  } else {
    done <- 0
    while (done == 0) {
      j <- which.max(lin_list$birth)
      # daughter <- lin_list$id[j] # nolint
      parent <- lin_list$parent[j]
      parent_j <- which(parent == lin_list$id)
      parent_in_list <- length(parent_j)
      if (parent_in_list == 1) {
        spec1 <- paste0(
          lin_list$label[parent_j],
          ":",
          lin_list$t_end[parent_j] - lin_list$birth[j]
        )
        spec2 <- paste0(
          lin_list$label[j],
          ":",
          lin_list$t_end[j] - lin_list$birth[j]
        )
        lin_list$label[parent_j] <- paste0("(", spec1, ",", spec2, ")")
        lin_list$t_end[parent_j] <- lin_list$birth[j]
        brts <- c(brts, lin_list$birth[j])
        lin_list <- lin_list[-j, ]
      } else {
        lin_list[j, 1:3] <- l_0_redux[which(l_0_redux[, 3] == parent), 1:3]
      }
      if (nrow(lin_list) == 1) {
        done <- 1
      }
    }
    brts <- rev(sort(age - brts))
    if (n_0 == 1) {
      brts <- c(age, brts)
    }
  }
  brts
}
