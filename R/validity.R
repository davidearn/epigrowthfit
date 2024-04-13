isTrueFalse <-
function(x)
	is.logical(x) && length(x) == 1L && !is.na(x)

isNumber <-
function(x)
	is.numeric(x) && length(x) == 1L && is.finite(x)

isString <-
function(x)
	is.character(x) && length(x) == 1L && !is.na(x)

isFlag <-
function(x) # not NA_integer_ when coerced to "integer"
	(is.logical(x) || is.numeric(x)) && length(x) == 1L && is.finite(x) &&
		(!is.double(x) || abs(x) < 0x1p+31)
