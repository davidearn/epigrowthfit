test_that("coef", {
  o <- egf_cache("egf-2.rds")
  co0 <- coef(o, full = FALSE)
  co1 <- coef(o, full = TRUE)

})

test_that("fixef", {

})

test_that("ranef", {

})

test_that("vcov", {

})

test_that("getCall", {
  o <- list(call = call("egf.method"))
  class(o) <- "egf"
  expect_identical(getCall(o), call("egf"))
  class(o) <- "egf_no_fit"
  expect_identical(getCall(o), call("egf"))
})

test_that("model.frame", {

})

test_that("model.matrix", {

})

test_that("terms", {

})

test_that("formula", {

})

test_that("nobs", {

})

test_that("df.residual", {

})

test_that("logLik", {

})

test_that("extractAIC", {

})
