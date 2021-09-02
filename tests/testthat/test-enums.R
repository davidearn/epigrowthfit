test_that("egf_get_flag", {
  flag_curve <- egf_get_flag("curve", c("exponential", "invalid name"))
  expect_type(flag_curve, "integer")
  expect_length(flag_curve, 2L)
  expect_gte(flag_curve[[1L]], 0L)
  expect_equal(flag_curve[[2L]], -1L)

  flag_family <- egf_get_flag("family", c("pois", "invalid name"))
  expect_type(flag_family, "integer")
  expect_length(flag_family, 2L)
  expect_gte(flag_family[[1L]], 0L)
  expect_equal(flag_family[[2L]], -1L)

  flag_prior <- egf_get_flag("prior", c("norm", "invalid name"))
  expect_type(flag_prior, "integer")
  expect_length(flag_prior, 2L)
  expect_gte(flag_prior[[1L]], 0L)
  expect_equal(flag_prior[[2L]], -1L)

  flag_test <- egf_get_flag("test", c("logspace_diff", "invalid name"))
  expect_type(flag_test, "integer")
  expect_length(flag_test, 2L)
  expect_gte(flag_test[[1L]], 0L)
  expect_equal(flag_test[[2L]], -1L)

  expect_error(egf_get_flag("invalid name", c("foo", "bar")))
})
