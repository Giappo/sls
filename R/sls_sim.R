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
  l_2 = sls::sls_sim.get_standard_l_2(),
  l_matrix_size = 1e4
) {

  # create the parameters
  pars <- sls_sim.get_pars(
    lambdas = lambdas,
    mus = mus,
    Ks = Ks
  )
  good_sim <- 0
  while (!good_sim) {

    # initialize data
    data <- sls_sim.initialize_data_new_clade(clade = 0, l_2 = l_2); clade <- 1;
    for (clade in l_2$clade_id) {

      # initialize data for the clade
      data <- sls_sim.initialize_data_new_clade(
        data = data,
        clade = clade,
        pars = pars,
        l_2 = l_2,
        l_matrix_size = l_matrix_size
      )
      while (data$t[[clade]] > 0) {

        # sample delta_n and delta_t
        deltas <- sls_sim.sample_deltas(
          data = data,
          clade = clade,
          pars = pars
        ); deltas

        # decide the event
        event <- sls_sim.decide_event(
          data = data,
          clade = clade,
          l_2 = l_2,
          deltas = deltas
        ); event

        # modify data accordingly
        output <- sls_sim.use_event(
          data = data,
          clade = clade,
          l_2 = l_2,
          event = event,
          deltas = deltas
        ); output
        data <- output
      }
    }

    # is the simulation in agreement with the conditioning?
    good_sim <- sls_sim.conditioning(
      data = data,
      l_2 = l_2,
      cond = cond
    ); good_sim
  }

  # retrieve branching times info
  brts <- sls_sim.get_brts(
    data = data,
    l_2 = l_2
  )

  return(
    list(
      l_tables = data$l_1,
      brts = brts
    )
  )
}
