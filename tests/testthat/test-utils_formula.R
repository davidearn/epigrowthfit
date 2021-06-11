test_that("`negate()` negates names", {
  x <- quote(x)
  minus_x <- call("-", x)
  expect_identical(negate(x), minus_x)
  expect_identical(negate(minus_x), x)
})

test_that("`negate()` negates calls", {
  x <- quote(x)
  log_x <- call("log", x)
  minus_log_x <- call("-", log_x)
  expect_identical(negate(log_x), minus_log_x)
  expect_identical(negate(minus_log_x), log_x)
})

test_that("`negate()` negates atomic scalars", {
  expect_identical(negate(1), call("-", 1))
  expect_identical(negate(call("-", 1)), 1)
  expect_identical(negate(-1), call("-", -1))
  expect_identical(negate(call("-", -1)), -1)
})

test_that("`unsplit_terms()` and `split_terms()` construct and deconstruct nested calls to `+` and `-` not involving `|`", {
  x <- quote(1 + a * b - b)
  l <- list(1, quote(a * b), quote(-b))
  expect_identical(split_terms(x), l)
  expect_identical(unsplit_terms(l), x)
})

test_that("`unsplit_terms()` and `split_terms()` construct and deconstruct nested calls to `+` and `-` involving `|`", {
  x <- quote(w + (x | f/g) + (y + z | h))
  x_no_paren <- call("+", call("+", quote(w), quote(x | f/g)), quote(y + z | h))
  l <- list(quote(w), quote(x | f/g), quote(y + z | h))
  expect_identical(split_terms(x), l)
  expect_identical(split_terms(x_no_paren), l)
  expect_identical(unsplit_terms(l), x_no_paren)
})

test_that("`split_effects()` separates random effects terms from mixed effects formulae", {
  x <- y ~ x + (1 | f) + (a + b | g)
  l <- list(fixed = y ~ x, random = list(quote(1 | f), quote(a + b | g)))
  expect_identical(split_effects(x), l)
})

test_that("`split_interaction()` deconstructs nested calls to `:`", {
  x <- quote(a:b:I(f:g):log(h))
  l <- list(quote(a), quote(b), quote(I(f:g)), quote(log(h)))
  expect_identical(split_interaction(x), l)
})

test_that("`gsub_bar_plus()` replaces `|` with `+` in random effects terms", {
  x1 <- ~x + (1 | f) + (a + b | g)
  x2 <- ~x
  ## Expected result does not use `(` explicitly
  x2[[2L]] <- call("+", call("+", x2[[2L]], quote(1 + f)), quote(a + b + g))
  expect_identical(gsub_bar_plus(x1), x2)
})

test_that("`simplify_terms()` matches `terms(simplify = TRUE)` with fixed effects formulae", {
  x1 <- ~0 + x * y - x
  x2 <- formula(terms(x1, simplify = TRUE))
  expect_identical(simplify_terms(x1), x2)
})

test_that("`simplify_terms()` behaves as expected with mixed effects formulae", {
  x0 <- ~0 + x * y - y
  x1 <- ~0 + x * y - y + (1 | f/g) + (a | f) + (0 + b | f:g)
  x2 <- formula(terms(x0, simplify = TRUE))
  ## Expected result does not use `(` explicitly
  x2[[2L]] <- call("+", call("+", x2[[2L]], quote(a | f)), quote(b - 1 | f:g))
  expect_identical(simplify_terms(x1), x2)
})
