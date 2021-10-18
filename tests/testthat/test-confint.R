co <- readRDS(system.file("exdata", "confint-egf.rds", package = "epigrowthfit", mustWork = TRUE))

test_that("basic", {
  fo <- readRDS(system.file("exdata", "fitted-egf.rds", package = "epigrowthfit", mustWork = TRUE))
  expect_type(co, "list")
  expect_s3_class(co, c("egf_confint", "data.frame"), exact = TRUE)
  expect_identical(c(co), c(confint(fo)))
})

skip_on_cran()

test_that("plot", {
  bars <- function() {
    op <- par(mar = c(4.5, 4, 2, 1), oma = c(0, 0, 0, 0))
    on.exit(par(op))
    plot(cio, type = "bars")
  }
  expect_doppelganger("plot-egf_confint-bars", fig = bars)

  boxes <- function() {
    op <- par(mar = c(0.2, 0, 0.2, 0), oma = c(4.5, 6, 2, 1), las = 1)
    on.exit(par(op))
    plot(cio, type = "boxes")
  }
  expect_doppelganger("plot-egf_confint-boxes", fig = boxes)
})
