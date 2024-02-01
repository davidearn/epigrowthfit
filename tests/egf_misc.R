attach(asNamespace("epigrowthfit"))
library(methods)
library(tools)
options(warn = 2L, error = if (interactive()) recover)

o.1 <- egf_cache("egf-1.rds")
o.2 <- egf_cache("egf-2.rds")


## egf_top #############################################################

x1 <- NULL
x2 <- egf_model()
x3 <- list(model = x2)
class(x3) <- "egf"

s1 <- egf_top(x1, link = FALSE)
s2 <- egf_top(x2, link = FALSE)
s3 <- egf_top(x3, link = FALSE)

stopifnot(exprs = {
	is.character(s1)
	length(s1) > 0L
	!anyNA(s1)
	is.null(names(s1))

	is.character(s2)
	length(s2) > 0L
	match(s2, s1, 0L) > 0L

	identical(s3, s2)
})


## egf_has_random ######################################################

e <- new.env()
e[["data"]] <- list()
o <- list(tmb_out = list(env = e))
class(o) <- "egf"

e[["data"]][["Z"]] <- matrix(0, 3L, 3L)
stopifnot( egf_has_random(o))

e[["data"]][["Z"]] <- matrix(0, 3L, 0L)
stopifnot(!egf_has_random(o))


## egf_has_converged ###################################################

o <- list(optimizer_out = list(convergence = 0L),
          value = 100,
          gradient = runif(10L, -0.5, 0.5),
          hessian = TRUE)
class(o) <- c("egf", "list") # "list" only for dispatch to within.list

stopifnot(exprs = {
	egf_has_converged(o)
	!egf_has_converged(within(o, optimizer_out[["convergence"]] <- 1L))
	!egf_has_converged(within(o, value <- NaN))
	!egf_has_converged(within(o, gradient[1L] <- 10))
	!egf_has_converged(within(o, hessian <- FALSE))
	is.na(egf_has_converged(within(o, hessian <- NA)))
})


## egf_(expand|condense)_par ###########################################

t.2 <- o.2[["tmb_out"]]
p.2 <- t.2[["env"]][["last.par.best"]]

p.2e <- egf_expand_par(t.2, p.2)
p.2ec <- egf_condense_par(t.2, p.2e)

len <- c(table(factor(names(p.2), levels = c("beta", "theta", "b"))))

stopifnot(exprs = {
	all.equal(p.2e, structure(c(p.2[1:4], theta = 0, p.2[-(1:4)]),
	                          lengths = len + c(0L, 1L, 0L)))
	all.equal(p.2ec, structure(p.2, lengths = len))
})


## egf_get_sdreport ####################################################

identical. <-
function(x, y, ...) {
	x[["env"]] <- y[["env"]] <- NULL
	identical(x, y, ...)
}

sd.1 <- o.1[["sdreport"]]
sd.1. <- egf_get_sdreport(o.1)
stopifnot(identical.(sd.1., sd.1))

o.1[["sdreport"]] <- NULL
assertWarning(sd.1. <- egf_get_sdreport(o.1))
stopifnot(identical.(sd.1., sd.1))


## egf_preprofile ######################################################
## FIXME: only testing a trivial case as no elements of 'beta' are mapped

l <- egf_preprofile(o.1, subset = 1:20, top = c("log(r)", "log(c0)"))

stopifnot(exprs = {
	is.list(l)
	length(l) == 2L
	identical(names(l), c("Y", "A"))
	identical(unname(l[["Y"]]),
	          new("dgCMatrix",
	              Dim = c(20L, 2L),
	              p = c(0L, 0L, 0L)))
	identical(unname(l[["A"]]),
	          new("dgCMatrix",
	              Dim = c(40L, 5L),
	              p = c(0L, 20L, 40L, 40L, 40L, 40L),
	              i = 0:39,
	              x = rep.int(1, 40L)))
})


## egf_cache ###########################################################

subdir <- R_user_dir("epigrowthfit", "cache")
lf1 <- list.files(subdir)
file <- "test.rds"
a <- 1
x <- egf_cache(file, a)
y <- egf_cache(file, a <- 2)
lf2 <- list.files(subdir)
z <- egf_cache(file, clear = TRUE)
lf3 <- list.files(subdir)

stopifnot(exprs = {
	identical(x, 1)
	identical(y, 1)
	identical(z, 0L)
	identical(a, 1)
	identical(setdiff(lf2, lf1), file)
	identical(lf3, lf1)
})
