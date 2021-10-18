so <- readRDS(system.file("exdata", "summary-egf.rds", package = "epigrowthfit", mustWork = TRUE))

test_that("basic", {
  o <- readRDS(system.file("exdata", "egf.rds", package = "epigrowthfit", mustWork = TRUE))
  fo <- readRDS(system.file("exdata", "fitted-egf.rds", package = "epigrowthfit", mustWork = TRUE))

  expect_type(so, "list")
  expect_s3_class(so, "egf_summary")
  expect_length(so, 5L)
  expect_named(so, c("fitted", "convergence", "value", "gradient", "hessian"), ignore.order = TRUE)
  expect_identical(so$convergence, o$optimizer_out$convergence)
  expect_identical(so$value, o$value)
  expect_identical(so$gradient, o$gradient)
  expect_identical(so$hessian, o$hessian)

  expect_type(so$fitted, "double")
  expect_identical(dim(so$fitted), c(6L, nlevels(fo$top)))
  expect_identical(dimnames(so$fitted), list(names(summary(0)), levels(fo$top)))
  for (s in colnames(so$fitted)) {
    eval(bquote(expect_equal(so$fitted[, .(s)], c(summary(fo$estimate[fo$top == .(s)])))))
  }
})

test_that("print", {
  capture.output({
    pso <- expect_condition(print(so), regexp = NA)
    expect_identical(pso, so)
    expect_invisible(print(so))
  })
})
