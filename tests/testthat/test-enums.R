test_that("egf_get_flag", {
  f <- function(type, enum) egf_get_flag(c(enum, "invalid enum"), type)
  flag <- Map(f,
    type = c("curve", "family", "prior"),
    enum = c("exponential", "pois", "norm")
  )
  for (s in names(flag)) {
    eval(bquote({
      expect_type(flag[[.(s)]], "integer")
      expect_length(flag[[.(s)]], 2L)
      expect_gte(flag[[.(s)]][1L], 0L)
      expect_identical(flag[[.(s)]][2L], -1L)
    }))
  }
  expect_error(egf_get_flag(c("foo", "bar"), "invalid type"))
})
