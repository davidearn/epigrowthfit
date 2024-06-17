library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)
example("egf", package = "epigrowthfit"); o.1 <- m1; o.2 <- m2

.S3method("all.equal", "confint.egf",
          function(target, current, ignore = NULL, ...) {
          	if (!is.null(ignore))
          		attributes(target)[ignore] <- attributes(current)[ignore] <-
          			NULL
          	NextMethod()
          })


## object ##############################################################

o.1c.w <- confint(o.1, A = NULL, method =    "wald", class = TRUE,
                  random = TRUE)
o.1c.p <- confint(o.1, A = NULL, method = "profile", class = TRUE,
                  top = "log(r)", subset = quote(country == "A" & wave == 1))
o.1c.u <- confint(o.1, A = NULL, method = "uniroot", class = TRUE,
                  top = "log(r)", subset = quote(country == "A" & wave == 1))

o.1f <- fitted(o.1, class = TRUE, se = TRUE)
o.1fc <- confint(o.1f, class = TRUE)

o.1p <- profile(o.1, A = NULL,
                top = "log(r)", subset = quote(country == "A" & wave == 1))
o.1pc <- confint(o.1p, class = TRUE)

stopifnot(exprs = {
	is.list(o.1c.w)
	identical(oldClass(o.1c.w), c("confint.egf", "data.frame"))
	all.equal(o.1c.w, o.1fc)

	is.list(o.1c.p)
	identical(oldClass(o.1c.p), c("confint.egf", "data.frame"))
	all.equal(o.1c.p, o.1pc)

	is.list(o.1c.u)
	identical(oldClass(o.1c.u), c("confint.egf", "data.frame"))
	all.equal(o.1c.u, o.1c.p, tolerance = 5e-06)
})


## parallel ############################################################

f <-
function(method, cores)
	confint(o.1, A = NULL, method = "uniroot", class = TRUE,
	        top = "log(r)", subset = quote(country == "A" & wave == 1),
	        parallel = egf_parallel(method = method, cores = cores))

windows <- .Platform[["OS.type"]] == "windows"
stopifnot(exprs = {
	all.equal(o.1c.u, f("multicore", if (windows) 1L else 2L))
	all.equal(o.1c.u, f("snow", 2L))
})


## plot ################################################################

op <- par(mar = c(4.5, 4, 2, 1), oma = c(0, 0, 0, 0))
plot(o.1c.w)
par(op)
