library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)
example("egf", package = "epigrowthfit"); o.1 <- m1; o.2 <- m2


## print ###############################################################

vv <- withVisible(print(o.1))
stopifnot(exprs = {
	identical(vv[["value"]], o.1)
	identical(vv[["visible"]], FALSE)
})
