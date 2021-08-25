test_that("compute_tdoubling", {
  r <- c(rlnorm(10L, -3, 1), 0, NA, NaN, Inf)
  per <- 1L
  tdoubling <- compute_tdoubling(r = r, per = per)

  expect_equal(as.numeric(tdoubling), log(2) / r)
  expect_s3_class(tdoubling, "tdoubling")
  expect_equal(attr(tdoubling, "per"), per)
  expect_warning(compute_tdoubling(-1), "NA")

  expect_equal(print(tdoubling), as.numeric(tdoubling))
  expect_invisible(print(tdoubling))
})
