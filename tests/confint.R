library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)

.S3method("all.equal", "egf_confint",
          function(target, current, ignore = NULL, ...) {
          	if (!is.null(ignore))
          		attributes(target)[ignore] <- attributes(current)[ignore] <-
          			NULL
          	NextMethod()
          })

o.1    <- egf_cache(        "egf-1.rds")
o.1c.w <- egf_cache("confint-egf-1.rds") # confint(method =    "wald")
o.1c.p <- egf_cache("confint-egf-2.rds") # confint(method = "profile")
o.1c.u <- egf_cache("confint-egf-3.rds") # confint(method = "uniroot")
o.1f   <- egf_cache( "fitted-egf-1.rds")
o.1p   <- egf_cache("profile-egf-1.rds")


## object ##############################################################

o.1fc <- confint(o.1f)
o.1pc <- confint(o.1p)
o.1pc[["linear_combination"]] <- NULL

stopifnot(exprs = {
	is.list(o.1c.w)
	identical(oldClass(o.1c.w), c("egf_confint", "data.frame"))
	all.equal(o.1c.w, o.1fc, ignore = c("class", "method", "frame_windows"))

	is.list(o.1c.p)
	identical(oldClass(o.1c.p), c("egf_confint", "data.frame"))
	all.equal(o.1c.p, o.1pc, ignore = c("class", "method", "frame_windows", "A", "x"))

	is.list(o.1c.u)
	identical(oldClass(o.1c.u), c("egf_confint", "data.frame"))
	all.equal(o.1c.u, o.1c.p, ignore = "method", tolerance = 1e-03)
})


## parallel ############################################################

f <-
function(method)
	confint(o.1,
	        method = "uniroot",
	        subset = country == "A" & wave == 1,
	        parallel = egf_parallel(method = method, cores = 2L))

stopifnot(exprs = {
	all.equal(o.1c.u, f("multicore"))
	all.equal(o.1c.u, f("snow"))
})


## plot ################################################################

op <- par(mar = c(4.5, 4, 2, 1), oma = c(0, 0, 0, 0))
plot(o.1c.p, type = "bars")
par(op)

op <- par(mar = c(0.2, 0, 0.2, 0), oma = c(4.5, 6, 2, 1), las = 1)
plot(o.1c.p, type = "boxes")
par(op)
