attach(asNamespace("epigrowthfit"))
library(tools)
options(warn = 2L, error = if (interactive()) recover)
example("egf", package = "epigrowthfit"); o.1 <- m1; o.2 <- m2


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


## egf_par_expand   ####################################################
## egf_par_condense ####################################################

t.2 <- o.2[["tmb_out"]]
p.2 <- t.2[["env"]][["last.par.best"]]

p.2e <- egf_par_expand(t.2, p.2)
p.2ec <- egf_par_condense(t.2, p.2e)

len <- c(table(factor(names(p.2), levels = c("beta", "theta", "b"))))

stopifnot(exprs = {
	all.equal(p.2e, structure(c(p.2[1:4], theta = 0, p.2[-(1:4)]),
	                          len = len + c(0L, 1L, 0L), names = NULL))
	all.equal(p.2ec, structure(p.2, len = len, names = NULL))
})


## egf_report   ########################################################
## egf_adreport ########################################################

identical. <-
function(x, y, ...) {
	x[["env"]] <- y[["env"]] <- NULL
	identical(x, y, ...)
}

r.1 <- o.1[["tmb_out"]][["env"]][[".__egf__."]][["report"]]
r.1. <- egf_report(o.1)
stopifnot(identical.(r.1., r.1))

o.1[["tmb_out"]][["env"]][[".__egf__."]][["report"]] <- NULL
r.1. <- egf_report(o.1)
stopifnot(identical.(r.1., r.1))

r.1 <- o.1[["tmb_out"]][["env"]][[".__egf__."]][["adreport"]]
r.1. <- egf_adreport(o.1)
stopifnot(identical.(r.1., r.1))

o.1[["tmb_out"]][["env"]][[".__egf__."]][["adreport"]] <- NULL
assertCondition(r.1. <- egf_adreport(o.1), "message")
stopifnot(identical.(r.1., r.1))
