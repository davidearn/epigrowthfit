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


## isFlag ##############################################################

stopifnot(exprs = {
	 isFlag(TRUE)
	 isFlag(FALSE)
	!isFlag(NA)
	 isFlag(1L)
	 isFlag(-1)
	 isFlag(-1.1)
	!isFlag(0x1p+31)
	!isFlag(NA_real_)
	!isFlag(NaN)
	!isFlag(Inf)
	!isFlag(NULL)
	!isFlag(c(1, 1))
	!isFlag("1")
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
