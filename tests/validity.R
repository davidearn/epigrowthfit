attach(asNamespace("epigrowthfit"))
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## isTrueFalse #########################################################

stopifnot(exprs = {
	 isTrueFalse(TRUE)
	 isTrueFalse(FALSE)
	!isTrueFalse(NA)
	!isTrueFalse(NULL)
	!isTrueFalse(c(TRUE, TRUE))
	!isTrueFalse(1)
})


## isInteger ###########################################################

stopifnot(exprs = {
	 isInteger(TRUE)
	 isInteger(FALSE)
	!isInteger(NA)
	 isInteger(1L)
	 isInteger(-1)
	 isInteger(-1.1)
	!isInteger(0x1p+31)
	!isInteger(NA_real_)
	!isInteger(NaN)
	!isInteger(Inf)
	!isInteger(NULL)
	!isInteger(c(1, 1))
	!isInteger("1")
})


## isNumber ############################################################

stopifnot(exprs = {
	!isNumber(TRUE)
	 isNumber(1L)
	 isNumber(-1)
	 isNumber(-1.1)
	 isNumber(0x1p+31)
	!isNumber(NA_real_)
	!isNumber(NaN)
	!isNumber(Inf)
	!isNumber(NULL)
	!isNumber(c(1, 1))
	!isNumber("1")
})


## isString ############################################################

stopifnot(exprs = {
	 isString("1")
	!isString(NA_character_)
	!isString(NULL)
	!isString(c("1", "1"))
	!isString(1)
})
