attach(asNamespace("epigrowthfit"))
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## egf_enum ############################################################

ii <- Map(function(enum, type) egf_enum(c(enum, "invalid enum"), type),
          enum = c("exponential", "pois", "norm"),
          type = c("curve", "family", "prior"))

for (i in ii)
	stopifnot(exprs = {
		is.integer(i)
		length(i) == 2L
		i[1L] >=  0L
		i[2L] == -1L
	})
assertError(egf_enum(c("foo", "bar"), "invalid type"))
