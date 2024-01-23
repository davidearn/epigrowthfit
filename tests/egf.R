library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)

o.1 <- egf_cache("egf-1.rds")


## print ###############################################################

vv <- withVisible(print(o.1))
stopifnot(exprs = {
	identical(vv[["value"]], o.1)
	identical(vv[["visible"]], FALSE)
})
