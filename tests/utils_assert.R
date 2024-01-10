library(epigrowthfit)
options(warn = 2L, error = recover)


## stopifnot1 ######
    expect_null(stopifnot1(TRUE))
    assertError(stopifnot1(FALSE))
    assertError(stopifnot1(NA))
    assertError(stopifnot1(NULL))
    assertError(stopifnot1(c(TRUE, FALSE)))
    assertError(stopifnot1(1))
    assertError(stopifnot1(TRUE, FALSE))
    assertError(stopifnot1(FALSE, m = "bar"), "bar")


## warnifnot1 ######
    expect_null(warnifnot1(TRUE))
    assertWarning(warnifnot1(FALSE))
    assertWarning(warnifnot1(NA))
    assertWarning(warnifnot1(NULL))
    assertWarning(warnifnot1(c(TRUE, FALSE)))
    assertWarning(warnifnot1(1))
    assertWarning(warnifnot1(TRUE, FALSE))
    assertWarning(warnifnot1(FALSE, m = "bar"), "bar")


## is_true_or_false ######
    expect_true(is_true_or_false(TRUE))
    expect_true(is_true_or_false(FALSE))
    expect_false(is_true_or_false(NA))
    expect_false(is_true_or_false(NULL))
    expect_false(is_true_or_false(c(TRUE, TRUE)))
    expect_false(is_true_or_false(1))


## is_flag ######
    expect_true(is_flag(TRUE))
    expect_true(is_flag(FALSE))
    expect_false(is_flag(NA))
    expect_true(is_flag(1))
    expect_true(is_flag(0))
    expect_true(is_flag(-1))
    expect_false(is_flag(NA_real_))
    expect_false(is_flag(NULL))
    expect_false(is_flag(c(1, 1)))


## is_number ######
    expect_true(is_number(1L))
    expect_true(is_number(1))
    expect_true(is_number(1e+10))
    expect_true(is_number(1.1))
    expect_false(is_number(NA_real_))
    expect_false(is_number(NaN))
    expect_false(is_number(Inf))
    expect_false(is_number(NULL))
    expect_false(is_number(c(1, 1)))
    expect_false(is_number(TRUE))
    expect_false(is_number("1"))

    expect_true(is_number(1, kind = "positive"))
    expect_false(is_number(0, kind = "positive"))
    expect_false(is_number(-1, kind = "positive"))

    expect_true(is_number(1, kind = "nonnegative"))
    expect_true(is_number(0, kind = "nonnegative"))
    expect_false(is_number(-1, kind = "nonnegative"))

    expect_false(is_number(1, kind = "negative"))
    expect_false(is_number(0, kind = "negative"))
    expect_true(is_number(-1, kind = "negative"))

    expect_false(is_number(1, kind = "nonpositive"))
    expect_true(is_number(0, kind = "nonpositive"))
    expect_true(is_number(-1, kind = "nonpositive"))

    expect_true(is_number(1L, integer = TRUE))
    expect_true(is_number(1, integer = TRUE))
    expect_false(is_number(1 + .Machine$double.eps, integer = TRUE))


## is_number_in_interval ######
    expect_false(is_number_in_interval(0, 0, 1, "()"))
    expect_false(is_number_in_interval(0, 0, 1, "(]"))
    expect_true(is_number_in_interval(0, 0, 1, "[)"))
    expect_true(is_number_in_interval(0, 0, 1, "[]"))



## is_string ######
    expect_true(is_string("foo"))
    expect_false(is_string(NA_character_))
    expect_false(is_string(NULL))
    expect_false(is_string(c("foo", "foo")))
    expect_false(is_string(1))

