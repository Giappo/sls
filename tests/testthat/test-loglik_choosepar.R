context("loglik_choosepar")

test_that("use", {

  skip("function not in use")

  brts_m <- c(3, 2, 1)
  brts_s <- c(2.5, 1.5)
  pars2 <- c(100, 1, 1, 1, 0, 2)
  idparsopt <- 4
  idparsfix <- 1:3
  trparsopt <- c(0.2, 0.1, 0.4, 0.2)[idparsopt]
  trparsfix <- c(0.2, 0.1, 0.4, 0.2)[idparsfix]
  sls_loglik_choosepar(
    trparsopt = trparsopt,
    trparsfix = trparsfix,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    idparsnoshift = NULL,
    pars2 = pars2,
    brts_m = brts_m,
    brts_s = brts_s
  )

})

test_that("abuse", {

  skip("function not in use")
})
