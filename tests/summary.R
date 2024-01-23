library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)

o.1  <- egf_cache(        "egf-1.rds")
o.1s <- egf_cache("summary-egf-1.rds")
o.1f <- egf_cache( "fitted-egf-1.rds")


## object ##############################################################

## FIXME: working around apparent bug in simplify2array
s2a <-
function(...) {
	r <- simplify2array(...)
	array(r, dim = unname(dim(r)), dimnames = dimnames(r))
}

stopifnot(exprs = {
	is.list(o.1s)
	identical(oldClass(o.1s), "egf_summary")
	length(o.1s) == 5L
	identical(names(o.1s), c("fitted", "convergence", "value", "gradient", "hessian"))
	identical(o.1s[["convergence"]], o.1[["optimizer_out"]][["convergence"]])
	identical(o.1s[["value"]], o.1[["value"]])
	identical(o.1s[["gradient"]], o.1[["gradient"]])
	identical(o.1s[["hessian"]], o.1[["hessian"]])
	all.equal(o.1s[["fitted"]],
	          s2a(tapply(o.1f[["estimate"]], o.1f[["top"]], summary)))
})


## print ###############################################################

vv <- withVisible(print(o.1s))
stopifnot(exprs = {
	identical(vv[["value"]], o.1s)
	identical(vv[["visible"]], FALSE)
})
