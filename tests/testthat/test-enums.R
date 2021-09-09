test_that("egf_get_flag", {
  f <- function(type, enum) egf_get_flag(type, c(enum, "invalid name"))
  flag <- Map(f,
    type = c("curve", "family", "prior"),
    enum = c("exponential", "pois", "norm")
  )
  for (s in names(flag)) {
    eval(bquote({
      expect_type(flag[[.(s)]], "integer")
      expect_length(flag[[.(s)]], 2L)
      expect_gte(flag[[.(s)]][1L], 0L)
      expect_equal(flag[[.(s)]][2L], -1L)
    }))
  }
  expect_error(egf_get_flag("invalid name", c("foo", "bar")))
})
