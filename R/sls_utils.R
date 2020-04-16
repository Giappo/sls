#---- package specific functions
#' @title Transition matrix builder
#' @author Giovanni Laudanno
#' @description Builds the transition matrix to integrate
#'  the differential equations of the P-equation
#' @inheritParams default_params_doc
#' @return The transition matrix
#' @export
p_transition_matrix <- function(
  lambda,
  mu,
  matrix_size
) {
  nvec <- 0:(matrix_size - 1)
  m <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  m[row(m) == col(m) + 1] <- m[row(m) == col(m) + 1] +
    lambda * (nvec[1:(matrix_size - 1)])
  m[row(m) == col(m) - 1] <- m[row(m) == col(m) - 1] +
    mu * (nvec[2:matrix_size])
  m[row(m) == col(m)] <- m[row(m) == col(m)] -
    (lambda + mu) * (nvec[1:matrix_size])
  m[matrix_size, matrix_size] <- -m[matrix_size - 1, matrix_size]; m
  testit::assert(colSums(m) < 1e-10)
  return(m)
}

#' @title Check input
#' @author Giovanni Laudanno
#' @description It checks the inputs
#' @inheritParams default_params_doc
#' @return nothing
#' @export
sls_check_input <- function(
  brts_m,
  brts_s,
  cond,
  n_0,
  n_max
) {
  if (length(brts_m) <= 0) {
    stop("main clade branching times cannot be an empty vector")
  }
  if (length(brts_s) <= 0) {
    stop("sub clade branching times cannot be an empty vector")
  }
  if (any(brts_m < 0)) {
    stop("all the branching times for the main clade have to be non negative")
  }
  if (any(brts_s < 0)) {
    stop("all the branching times for the sub clade have to be non negative")
  }
  if (!(cond %in% sls_conds())) {
   stop("this conditioning is not implemented")
  }
  if (!(n_0 %in% sls_n_0s())) {
   stop("this n_0 is not implemented")
  }
  if (n_max <= 0) {
   stop("it's not going to work with maximum species set to 0 or less")
  }
}

#' @title Conditionings
#' @author Giovanni Laudanno
#' @description Gives the conditionings accepted by sls
#' @inheritParams default_params_doc
#' @return the conditionings
#' @export
sls_conds <- function() {
  conds <- c(0, 1, 2, 3)
  conds
}

#' @title Starting species
#' @author Giovanni Laudanno
#' @description Gives the amount of starting species accepted by sls
#' @inheritParams default_params_doc
#' @return the possible n_0s
#' @export
sls_n_0s <- function() {
  n_0s <- c(2)
  n_0s
}

#' @title Logliks with division
#' @author Giovanni Laudanno
#' @description Get all the loglik functions with division
#' @inheritParams default_params_doc
#' @return loglik functions with division in sls
#' @export
sls_logliks_div <- function() {
  fun_list <- ls(paste0("package:", get_pkg_name())) # nolint internal function
  div_funs <- fun_list[sapply(
    fun_list, function(x)
      any(grepl("loglik_sls", x)) &
      !any(grepl("nodiv", x))
  )]
  div_funs
}

#' @title Logliks with no division
#' @author Giovanni Laudanno
#' @description Get all the loglik functions with no division
#' @inheritParams default_params_doc
#' @return loglik functions with no division in sls
#' @export
sls_logliks_nodiv <- function() {
  fun_list <- ls(paste0("package:", get_pkg_name())) # nolint internal function
  nodiv_funs <- fun_list[sapply(
    fun_list, function(x)
      any(grepl("loglik", x)) &
      any(grepl("nodiv", x)) &
      !any(grepl("logliks", x))
  )]
  nodiv_funs
}

#' @title Logliks for the experiment
#' @author Giovanni Laudanno
#' @description Get the loglik functions to use for the experiment
#' @inheritParams default_params_doc
#' @return loglik functions to use for the experiment
#' @export
sls_logliks_experiment <- function() {
  fun_list <- ls(paste0("package:", get_pkg_name())) # nolint internal function
  sls_p_funs <- fun_list[sapply(
    fun_list, function(x)
      any(grepl("sls_p", x)) &
      !any(grepl("2", x)) &
      !any(grepl("dummy", x))
  )]
  sls_p_funs
}

#' @title Get package name
#' @author Giovanni Laudanno
#' @description Get package name
#' @inheritParams default_params_doc
#' @return Package name
#' @export
get_pkg_name <- function() {
  pkg_name <- "sls"
  pkg_name
}

#' Get the names of the parameters used in the sls model
#' @author Giovanni Laudanno
#' @export
get_param_names <- function() {
  c("lambda_m", "mu_m", "lambda_s", "mu_s")
}

#---- general functions
#' @title cat2
#' @author Giovanni Laudanno
#' @description If verbose == TRUE cats the message, otherwise stays silent
#' @inheritParams default_params_doc
#' @return prints on screen
#' @export
cat2 <- function(
  message,
  verbose
) {
 if (verbose == TRUE) {
   cat(message)
 } else {
   return()
 }
}

#' @title Transform parameters
#' @description Transform parameters according to y = x / (1 + x)
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @return transformed parameters
#' @export
pars_transform_forward <- function(pars) {
  pars <- as.numeric(unlist(pars))
  pars_transformed <- pars / (1 + pars)
  pars_transformed[which(pars == Inf)] <- 1
  pars_transformed
}

#' @title Transform parameters back
#' @description Transform parameters back according to x = y / (1 + y)
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @return the original parameters
#' @export
pars_transform_back <- function(pars_transformed) {
  pars_transformed <- as.numeric(unlist(pars_transformed))
  pars <- pars_transformed / (1 - pars_transformed)
  pars
}

#' @title Cut word "loglik" from a name
#' @author Giovanni Laudanno
#' @description Cut word "loglik" from a name
#' @inheritParams default_params_doc
#' @return clean name
#' @export
cut_loglik_from_name <- function(
  function_name
) {
  if (grepl("loglik", function_name)) {
    model_name <- gsub(
      "_loglik",
      "",
      gsub(
        "loglik_",
        "",
        function_name
      )
    )
  } else {
    model_name <- NA
  }
  model_name
}

#' @title Get function names
#' @author Giovanni Laudanno
#' @description Get function names
#' @inheritParams default_params_doc
#' @return function names
#' @export
get_function_names <- function(
  loglik_functions
) {
  pkg_name <- get_pkg_name() # nolint internal function
  fun_list <- ls(paste0("package:", pkg_name))
  error_message <- paste0(
    "This is not a likelihood function provided by ",
    pkg_name,
    "!"
  )

  if (!is.list(loglik_functions)) {
    if (length(loglik_functions) == 1) {
      loglik_functions <- list(loglik_functions)
    }
  }

  function_names <- rep(NA, length(loglik_functions))
  model_names <- which_function <- function_names
  for (m in seq_along(loglik_functions)) {
    if (is.character(loglik_functions[[m]])) {
      if (length(
        (find_function <- which(fun_list == loglik_functions[[m]]))
      ) == 0) {
        stop(error_message)
      }
      which_function[m] <- find_function
    } else {
      fun <- eval(loglik_functions[m])[[1]]
      for (i in seq_along(fun_list)) {
        if (all.equal(get(fun_list[i]), fun) == TRUE) {
          which_function[m] <- i
        }
      }
    }
    if (is.null(which_function[m]) | is.na(which_function[m])) {
      stop(error_message)
    }
    function_names[m] <- toString(fun_list[which_function[m]])
    model_names[m] <- cut_loglik_from_name(function_names[m])
  }
  function_names <- fun_list[which_function]
  model_names <- unname(
    sapply(function_names, FUN = cut_loglik_from_name) # nolint internal function
  )
  if (any(is.na(model_names))) {
    stop(error_message)
  }
  invisible(function_names)
}

#' @title Check if provided models make sense
#' @author Giovanni Laudanno
#' @description Check if provided models make sense
#' @inheritParams default_params_doc
#' @return models names
#' @export
get_model_names <- function(
  function_names,
  verbose = FALSE
) {
  model_names <- function_names
  for (m in seq_along(function_names)) {
    model_names[m] <- cut_loglik_from_name(function_names[m])
    if (is.null(model_names[m]) | is.na(model_names[m])) {
      stop(paste0(
        "This is not a likelihood function provided by ",
        get_pkg_name(), # nolint internal function
        "!"
      ))
    }
  }
  if (verbose == TRUE) {
    cat("You are using the functions:", model_names)
  }
  model_names
}

#' @title Builds the right hand side of the ODE set
#' @description Builds the right hand side of the ODE set
#' @param time time
#' @param x vector to integrate
#' @param params transition matrix
#' @details This is not to be called by the user.
#' @author Giovanni Laudanno
#' @export
sls_loglik_rhs <- function(t, x, params) {
  list(params %*% x)
}
