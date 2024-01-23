library(epigrowthfit)
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## egf_get_flag ######
f <- function(type, enum) egf_get_flag(c(enum, "invalid enum"), type)
flag <- Map(f,
            type = c("curve", "family", "prior"),
            enum = c("exponential", "pois", "norm"))
for (s in names(flag)) {
    eval(bquote({
        is.integer(flag[[.(s)]])
        length(flag[[.(s)]]) == 2L
        flag[[.(s)]][1L] >= 0L
        identical(flag[[.(s)]][2L], -1L)
    }))
}
assertError(egf_get_flag(c("foo", "bar"), "invalid type"))

