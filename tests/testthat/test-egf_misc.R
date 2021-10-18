o <- readRDS(system.file("exdata", "egf.rds", package = "epigrowthfit", mustWork = TRUE))

test_that("egf_get_names_top", {
  names_top_all <- egf_get_names_top(NULL, link = FALSE)
  expect_type(names_top_all, "character")
  expect_gt(length(names_top_all), 0L)
  expect_false(anyNA(names_top_all))
  expect_named(names_top_all, NULL)

  model <- egf_model()
  names_top <- egf_get_names_top(model, link = FALSE)
  expect_type(names_top, "character")
  expect_gt(length(names_top), 0L)
  expect_true(all(names_top %in% names_top_all))

  object <- list(model = model)
  class(object) <- "egf"
  names_top_again <- egf_get_names_top(object, link = FALSE)
  expect_identical(names_top_again, names_top)
})

test_that("egf_has_random", {
  e <- new.env()
  e$data <- list()
  object <- list(tmb_out = list(env = e))
  class(object) <- "egf"

  e$data$Z <- matrix(numeric(9L), 3L, 3L)
  expect_true(egf_has_random(object))

  e$data$Z <- matrix(numeric(0L), 3L, 0L)
  expect_false(egf_has_random(object))
})

test_that("egf_has_converged", {
  object <- list(
    optimizer_out = list(convergence = 0L),
    value = 100,
    gradient = runif(10, -0.5, 0.5),
    hessian = TRUE
  )
  class(object) <- c("egf", "list")

  expect_true(egf_has_converged(object))
  expect_false(egf_has_converged(within(object, optimizer_out$convergence <- 1L)))
  expect_false(egf_has_converged(within(object, value <- NaN)))
  expect_false(egf_has_converged(within(object, gradient[1L] <- 10)))
  expect_false(egf_has_converged(within(object, hessian <- FALSE)))
  expect_identical(egf_has_converged(within(object, hessian <- NA)), NA)
})

test_that("egf_(expand|condense)_par", {
  ## FIXME: Trivial test because in example no parameters are mapped
  par <- o$tmb_out$env$last.par.best
  len <- c(table(factor(names(par), levels = c("beta", "theta", "b"))))
  epar <- egf_expand_par(o$tmb_out, par = par)
  expect_identical(epar, structure(par, lengths = len))
  cpar <- egf_condense_par(o$tmb_out, par = par)
  expect_identical(cpar, structure(par, lengths = len))
})

test_that("egf_get_sdreport", {
  sdr <- o$sdreport
  expect_identical(egf_get_sdreport(o)$cov.fixed, sdr$cov.fixed)
  o$sdreport <- NULL
  expect_warning(expect_equal(egf_get_sdreport(o)$cov.fixed, sdr$cov.fixed))
})

test_that("egf_preprofile", {
  ## FIXME: Trivial test because in example no parameters are mapped
  l <- egf_preprofile(o, subset = seq_len(20L), top = c("log(r)", "log(c0)"))
  expect_type(l, "list")
  expect_length(l, 2L)
  expect_named(l, c("Y", "A"), ignore.order = TRUE)
  expect_identical(unname(l$Y), Matrix(0, 20L, 2L))
  expect_identical(unname(l$A), sparseMatrix(i = seq_len(40L), j = rep.int(1:2, c(20L, 20L)), x = 1, dims = c(40L, 5L)))
})
