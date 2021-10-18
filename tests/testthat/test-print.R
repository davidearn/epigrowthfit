o <- readRDS(system.file("exdata", "egf.rds", package = "epigrowthfit", mustWork = TRUE))

test_that("basic", {
  capture.output({
    po <- expect_condition(print(o), regexp = NA)
    expect_identical(po, o)
    expect_invisible(print(o))
  })
})
