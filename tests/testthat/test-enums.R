test_that("get_flag", {
  flag_curve <- get_flag("curve", c("exponential", "invalid name"))
  expect_type(flag_curve, "integer")
  expect_length(flag_curve, 2L)
  expect_equal(flag_curve[[2L]], NA_integer_)

  flag_family <- get_flag("family", c("pois", "invalid name"))
  expect_type(flag_curve, "integer")
  expect_length(flag_curve, 2L)
  expect_equal(flag_curve[[2L]], NA_integer_)

  flag_prior <- get_flag("prior", c("norm", "invalid name"))
  expect_type(flag_prior, "integer")
  expect_length(flag_prior, 2L)
  expect_equal(flag_prior[[2L]], NA_integer_)

  expect_error(get_flag("invalid name", c("foo", "bar")))
})
