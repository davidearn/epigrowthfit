test_that("stop_if_not", {
  expect_null(stop_if_not(TRUE))
  expect_error(stop_if_not(FALSE))
  expect_error(stop_if_not(NA))
  expect_error(stop_if_not(NULL))
  expect_error(stop_if_not(c(TRUE, FALSE)))
  expect_error(stop_if_not(1))
  expect_error(stop_if_not(TRUE, FALSE))
  expect_error(stop_if_not(FALSE, m = "bar"), "bar")
})

test_that("warn_if_not", {
  expect_null(warn_if_not(TRUE))
  expect_warning(warn_if_not(FALSE))
  expect_warning(warn_if_not(NA))
  expect_warning(warn_if_not(NULL))
  expect_warning(warn_if_not(c(TRUE, FALSE)))
  expect_warning(warn_if_not(1))
  expect_warning(warn_if_not(TRUE, FALSE))
  expect_warning(warn_if_not(FALSE, m = "bar"), "bar")
})

test_that("stop_if_not_true_false", {
  expect_null(stop_if_not_true_false(TRUE))
  expect_null(stop_if_not_true_false(FALSE))
  expect_error(stop_if_not_true_false(NA))
  expect_error(stop_if_not_true_false(NULL))
  expect_error(stop_if_not_true_false(c(TRUE, TRUE)))
  expect_error(stop_if_not_true_false(1))

  expect_null(stop_if_not_true_false(1, allow_numeric = TRUE))
  expect_null(stop_if_not_true_false(0, allow_numeric = TRUE))
  expect_null(stop_if_not_true_false(-1, allow_numeric = TRUE))
  expect_error(stop_if_not_true_false(NA_real_, allow_numeric = TRUE))
  expect_error(stop_if_not_true_false(NULL, allow_numeric = TRUE))
  expect_error(stop_if_not_true_false(c(1, 1), allow_numeric = TRUE))
})

test_that("stop_if_not_integer", {
  expect_null(stop_if_not_integer(1L))
  expect_null(stop_if_not_integer(1))
  expect_null(stop_if_not_integer(1e+10))
  expect_error(stop_if_not_integer(1.1))
  expect_error(stop_if_not_integer(NA_integer_))
  expect_error(stop_if_not_integer(NaN))
  expect_error(stop_if_not_integer(Inf))
  expect_error(stop_if_not_integer(NULL))
  expect_error(stop_if_not_integer(c(1L, 1L)))
  expect_error(stop_if_not_integer(TRUE))
  expect_error(stop_if_not_integer("1"))

  expect_null(stop_if_not_integer(1L, kind = "positive"))
  expect_error(stop_if_not_integer(0L, kind = "positive"))
  expect_error(stop_if_not_integer(-1L, kind = "positive"))

  expect_null(stop_if_not_integer(1L, kind = "nonnegative"))
  expect_null(stop_if_not_integer(0L, kind = "nonnegative"))
  expect_error(stop_if_not_integer(-1L, kind = "nonnegative"))

  expect_error(stop_if_not_integer(1L, kind = "negative"))
  expect_error(stop_if_not_integer(0L, kind = "negative"))
  expect_null(stop_if_not_integer(-1L, kind = "negative"))

  expect_error(stop_if_not_integer(1L, kind = "nonpositive"))
  expect_null(stop_if_not_integer(0L, kind = "nonpositive"))
  expect_null(stop_if_not_integer(-1L, kind = "nonpositive"))
})

test_that("stop_if_not_number", {
  expect_null(stop_if_not_number(1L))
  expect_null(stop_if_not_number(1))
  expect_null(stop_if_not_number(1e+10))
  expect_null(stop_if_not_number(1.1))
  expect_error(stop_if_not_number(NA_real_))
  expect_error(stop_if_not_number(NaN))
  expect_error(stop_if_not_number(Inf))
  expect_error(stop_if_not_number(NULL))
  expect_error(stop_if_not_number(c(1, 1)))
  expect_error(stop_if_not_number(TRUE))
  expect_error(stop_if_not_number("1"))

  expect_null(stop_if_not_number(1, kind = "positive"))
  expect_error(stop_if_not_number(0, kind = "positive"))
  expect_error(stop_if_not_number(-1, kind = "positive"))

  expect_null(stop_if_not_number(1, kind = "nonnegative"))
  expect_null(stop_if_not_number(0, kind = "nonnegative"))
  expect_error(stop_if_not_number(-1, kind = "nonnegative"))

  expect_error(stop_if_not_number(1, kind = "negative"))
  expect_error(stop_if_not_number(0, kind = "negative"))
  expect_null(stop_if_not_number(-1, kind = "negative"))

  expect_error(stop_if_not_number(1, kind = "nonpositive"))
  expect_null(stop_if_not_number(0, kind = "nonpositive"))
  expect_null(stop_if_not_number(-1, kind = "nonpositive"))
})

test_that("stop_if_not_number_in_interval", {
  expect_error(stop_if_not_number_in_interval(0, 0, 1, "()"))
  expect_error(stop_if_not_number_in_interval(0, 0, 1, "(]"))
  expect_null(stop_if_not_number_in_interval(0, 0, 1, "[)"))
  expect_null(stop_if_not_number_in_interval(0, 0, 1, "[]"))
})

test_that("stop_if_not_string", {
  expect_null(stop_if_not_string("foo"))
  expect_error(stop_if_not_string(NA_character_))
  expect_error(stop_if_not_string(NULL))
  expect_error(stop_if_not_string(c("foo", "foo")))
  expect_error(stop_if_not_string(1))
})
