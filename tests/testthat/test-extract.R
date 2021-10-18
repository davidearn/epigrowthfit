test_that("getCall", {
  o <- list(call = call("egf.method"))
  class(o) <- "egf"
  expect_identical(getCall(o), call("egf"))
  class(o) <- "egf_no_fit"
  expect_identical(getCall(o), call("egf"))
})
