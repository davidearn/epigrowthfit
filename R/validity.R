isTrueFalse <-
function(x)
	is.logical(x) && length(x) == 1L && !is.na(x)

isInteger <- # TRUE if as.integer(x) has length 1, is not NA_integer_
function(x)
	(is.logical(x) || is.numeric(x)) && length(x) == 1L && is.finite(x) &&
		(!is.double(x) || abs(x) < 0x1p+31)

isNumber <-
function(x)
	is.numeric(x) && length(x) == 1L && is.finite(x)

isString <-
function(x)
	is.character(x) && length(x) == 1L && !is.na(x)
