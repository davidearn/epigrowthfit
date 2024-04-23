library(epigrowthfit)
library(methods)
options(warn = 2L, error = if (interactive()) recover)

o.1  <- egf_cache(        "egf-1.rds")
o.1p <- egf_cache("profile-egf-1.rds")


## object ##############################################################

stopifnot(exprs = {
	is.list(o.1p)
	identical(oldClass(o.1p), c("profile.egf", "profile"))
	length(o.1p) == 1L
	identical(names(o.1p), "log(r), A, window_01")
	identical(dim(o.1p), c(1L, 1L))
	identical(dimnames(o.1p), list("A, window_01", "log(r)"))

	is.list(o.1p[[1L]])
	identical(oldClass(o.1p[[1L]]), "data.frame")
	length(o.1p[[1L]]) == 2L
	identical(names(o.1p[[1L]]), c("z", "par.vals"))

	is.double(z <- o.1p[[1L]][["z"]])
	!is.matrix(z)
	min(abs(z)) == 0
	prod(sign(range(z))) == -1

	is.double(par.vals <- o.1p[[1L]][["par.vals"]])
	is.matrix(par.vals)
	ncol(par.vals) == 1L # for now
	!is.unsorted(par.vals, strictly = TRUE)
	par.vals[which.min(abs(z))] == coef(o.1)[1L]

	is.factor    (attr(o.1p, "top"   ))
	is.factor    (attr(o.1p, "ts"    ))
	is.factor    (attr(o.1p, "window"))
	is.data.frame(attr(o.1p, "frame" ))
	is           (attr(o.1p, "A"     ), "dgCMatrix")
	is.double    (attr(o.1p, "par"   ))
	identical    (attr(o.1p, "level" ), 0.95)
})


## confint #############################################################

o.1pc <- confint(o.1p, level = 0.95, class = TRUE)
n <- length(o.1p)

stopifnot(exprs = {
	is.list(o.1pc)
	identical(oldClass(o.1pc), c("confint.egf", "data.frame"))
	length(o.1pc) == 5L
	identical(names(o.1pc), c("top", "ts", "window", "value", "ci"))

	all(vapply(o.1pc[c("top", "ts", "window")], is.factor, FALSE))
	all(vapply(o.1pc[c("value", "ci"        )], is.double, FALSE))

	is.vector(o.1pc[["value"]])
	is.matrix(o.1pc[["ci"]])
	identical(dim(o.1pc[["ci"]]), c(1L, 2L))
	identical(dimnames(o.1pc[["ci"]]), list(NULL, c("2.5 %", "97.5 %")))
	all(o.1pc[["ci"]][, 1L] < o.1pc[["value"]])
	all(o.1pc[["ci"]][, 2L] > o.1pc[["value"]])
})


## parallel ############################################################

f <-
function(method, cores)
	profile(o.1, A = NULL,
	        top = "log(r)", subset = quote(country == "A" & wave == 1),
	        parallel = egf_parallel(method = method, cores = cores))

windows <- .Platform[["OS.type"]] == "windows"
stopifnot(exprs = {
	all.equal(o.1p, f("multicore", if (windows) 1L else 2L))
	all.equal(o.1p, f("snow", 2L))
})


## plot ################################################################

plot(o.1p, type = "z^2", bty = "u", las = 1)
