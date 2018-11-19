#' @title Simulate an sls process
#' @description Simulate an sls process
#' @inheritParams default_params_doc
#' @return  best parameters
#' @export
sls_sim <- function(
  lambdas,
  mus,
  Ks = c(Inf, Inf),
  cond = 3,
  LS = sls::sls_sim.get_standard_LS(),
  l_matrix_size = 1e4
) {
  pars <- sls_sim.get_pars(
    lambdas = lambdas,
    mus = mus,
    Ks = Ks
  )
  good_sim <- 0
  while (!good_sim) {
    data <- sls_sim.initialize_LL_new_clade(clade = 0); clade <- 1;
    for (clade in LS$clade_id) {
      data <- sls_sim.initialize_LL_new_clade(
        data = data,
        clade = clade,
        pars = pars,
        LS = LS,
        l_matrix_size = l_matrix_size
      )
      t <- sls_sim.initialize_t_new_clade(
        data = data,
        clade = clade
      )
      while (t > 0) {
        deltas <- sls_sim.sample_deltas(
          data = data,
          clade = clade,
          pars = pars
        ); delta_n <- deltas$delta_n; delta_t <- deltas$delta_t; deltas
        event <- sls_sim.decide_event(
          data = data,
          clade = clade,
          LS = LS,
          delta_n = delta_n,
          delta_t = delta_t,
          t = t
        )
        output <- sls_sim.use_event(
          data = data,
          clade = clade,
          LS = LS,
          event = event,
          t = t - deltas$delta_t
        )
        t <- output$t
        data <- output$data
      }
    }
    good_sim <- sls_sim.conditioning(
      data = data,
      LS = LS,
      cond = cond
    ); good_sim
  }

  brts <- sls_sim.get_brts(
    data = data,
    LS = LS
  )

  return(
    list(
      l_tables = data$LL,
      brts = brts
    )
  )
}
