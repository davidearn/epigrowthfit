attach(asNamespace("epigrowthfit"))
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## egf_get_flag ########################################################

flag <- Map(function(type, enum) egf_get_flag(c(enum, "invalid enum"), type),
            type = c("curve", "family", "prior"),
            enum = c("exponential", "pois", "norm"))

for (i in flag)
	stopifnot(exprs = {
		is.integer(i)
		length(i) == 2L
		i[1L] >= 0L
		i[2L] == -1L
	})
assertError(egf_get_flag(c("foo", "bar"), "invalid type"))
