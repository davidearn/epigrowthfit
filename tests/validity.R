attach(asNamespace("epigrowthfit"))
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## is_true_or_false ####################################################

stopifnot(exprs = {
	 is_true_or_false(TRUE)
	 is_true_or_false(FALSE)
	!is_true_or_false(NA)
	!is_true_or_false(NULL)
	!is_true_or_false(c(TRUE, TRUE))
	!is_true_or_false(1)
})


## is_flag #############################################################

stopifnot(exprs = {
	 is_flag(TRUE)
	 is_flag(FALSE)
	!is_flag(NA)
	 is_flag( 1)
	 is_flag( 0)
	 is_flag(-1)
	!is_flag(NA_real_)
	!is_flag(NULL)
	!is_flag(c(1, 1))
})


## is_number ###########################################################

stopifnot(exprs = {
	 is_number(1L)
	 is_number(1)
	 is_number(1e+10)
	 is_number(1.1)
	!is_number(NA_real_)
	!is_number(NaN)
	!is_number(Inf)
	!is_number(NULL)
	!is_number(c(1, 1))
	!is_number(TRUE)
	!is_number("1")

	 is_number( 1, kind = "positive")
	!is_number( 0, kind = "positive")
	!is_number(-1, kind = "positive")

	 is_number( 1, kind = "nonnegative")
	 is_number( 0, kind = "nonnegative")
	!is_number(-1, kind = "nonnegative")

	!is_number( 1, kind = "negative")
	!is_number( 0, kind = "negative")
	 is_number(-1, kind = "negative")

	!is_number( 1, kind = "nonpositive")
	 is_number( 0, kind = "nonpositive")
	 is_number(-1, kind = "nonpositive")

	 is_number(1L, integer = TRUE)
	 is_number(1, integer = TRUE)
	!is_number(1 + .Machine[["double.eps"]], integer = TRUE)
})


## is_number_in_interval ###############################################

stopifnot(exprs = {
	!is_number_in_interval(0, 0, 1, "()")
	!is_number_in_interval(0, 0, 1, "(]")
	 is_number_in_interval(0, 0, 1, "[)")
	 is_number_in_interval(0, 0, 1, "[]")
})


## is_string ###########################################################

stopifnot(exprs = {
	 is_string("foo")
	!is_string(NA_character_)
	!is_string(NULL)
	!is_string(c("foo", "foo"))
	!is_string(1)
})
