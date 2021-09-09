test_that("negate", {
  x <- quote(x)
  log_x <- call("log", x)

  minus_x <- call("-", x)
  expect_identical(negate(x), minus_x)
  expect_identical(negate(minus_x), x)

  minus_log_x <- call("-", log_x)
  expect_identical(negate(log_x), minus_log_x)
  expect_identical(negate(minus_log_x), log_x)

  expect_identical(negate(1), call("-", 1))
  expect_identical(negate(call("-", 1)), 1)
  expect_identical(negate(-1), call("-", -1))
  expect_identical(negate(call("-", -1)), -1)
})

test_that("(un)?split_terms", {
  x <- quote(1 + a * b - b)
  l <- list(1, quote(a * b), quote(-b))
  expect_identical(split_terms(x), l)
  expect_identical(unsplit_terms(l), x)

  x <- quote(w + (x | f/g) + (y + z | h))
  l <- list(quote(w), quote(x | f/g), quote(y + z | h))
  expect_identical(split_terms(x), l)
  expect_identical(unsplit_terms(l), x)
})

test_that("split_effects", {
  x <- y ~ x + (1 | f) + (a + b | g)
  l <- list(fixed = y ~ x, random = list(quote(1 | f), quote(a + b | g)))
  expect_identical(split_effects(x), l)
})

test_that("split_interaction", {
  x <- quote(a:b:I(f:g):log(h))
  l <- list(quote(a), quote(b), quote(I(f:g)), quote(log(h)))
  expect_identical(split_interaction(x), l)
})

test_that("gsub_bar_plus", {
  x1 <- ~x + (1 | f) + (a + b | g)
  x2 <- ~x
  ## Expected result does not use `(` explicitly
  x2[[2L]] <- call("+", call("+", x2[[2L]], quote(1 + f)), quote(a + b + g))
  expect_identical(gsub_bar_plus(x1), x2)
})

test_that("simplify_terms", {
  x1 <- ~0 + x * y - y
  y1 <- formula(terms(x1, simplify = TRUE))
  expect_identical(simplify_terms(x1), y1)

  x2 <- ~0 + x * y - y + (1 | f/g) + (a | f) + (0 + b | f:g)
  y2 <- y1
  y2[[2L]] <- call("+", call("+", y2[[2L]], quote((a | f))), quote((b - 1 | f:g)))
  expect_identical(simplify_terms(x2), y2)
})
