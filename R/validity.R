is_true_or_false <-
function(x)
	is.logical(x) && length(x) == 1L && !is.na(x)

is_flag <-
function(x)
	(is.logical(x) || is.numeric(x)) && length(x) == 1L && is.finite(x)

is_number <-
function(x,
         kind = c("any", "positive", "nonnegative", "negative", "nonpositive"),
         integer = FALSE) {
	kind <- match.arg(kind)
	relop <- switch(kind,
	                "any"         = function(e1, e2) TRUE,
	                "positive"    = `>`,
	                "nonnegative" = `>=`,
	                "negative"    = `<`,
	                "nonpositive" = `<=`)
	is.numeric(x) && length(x) == 1L && is.finite(x) &&
		relop(x, 0) && (!integer || x %% 1 == 0)
}

is_number_in_interval <-
function(x, a = -Inf, b = Inf, include = c("()", "(]", "[)", "[]")) {
	include <- match.arg(include)
	relop1 <- switch(substr(include, 1L, 1L), "(" = `>`, "[" = `>=`)
	relop2 <- switch(substr(include, 2L, 2L), ")" = `<`, "]" = `<=`)
	is.numeric(x) && length(x) == 1L && !is.na(x) &&
		relop1(x, a) && relop2(x, b)
}

is_string <-
function(x)
	is.character(x) && length(x) == 1L && !is.na(x)
