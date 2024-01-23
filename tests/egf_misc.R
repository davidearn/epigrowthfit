attach(asNamespace("epigrowthfit"))
library(methods)
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## egf_get_names_top ######
x <- egf_get_names_top(NULL, link = FALSE)
is.character(x)
length(x) > 0L
!anyNA(x)
is.null(names(x))

model <- egf_model()
x0 <- egf_get_names_top(model, link = FALSE)
is.character(x0)
length(x0) > 0L
all(x0 %in% x)

o <- list(model = model)
class(o) <- "egf"
identical(egf_get_names_top(o, link = FALSE), x0)


## egf_has_random ######
e <- new.env()
e$data <- list()
o <- list(tmb_out = list(env = e))
class(o) <- "egf"

e$data$Z <- matrix(double(9L), 3L, 3L)
egf_has_random(o)

e$data$Z <- matrix(double(0L), 3L, 0L)
!egf_has_random(o)


## egf_has_converged ######
o <- list(optimizer_out = list(convergence = 0L),
          value = 100,
          gradient = runif(10, -0.5, 0.5),
          hessian = TRUE)
class(o) <- c("egf", "list")

egf_has_converged(o)
!egf_has_converged(within(o, optimizer_out$convergence <- 1L))
!egf_has_converged(within(o, value <- NaN))
!egf_has_converged(within(o, gradient[1L] <- 10))
!egf_has_converged(within(o, hessian <- FALSE))
identical(egf_has_converged(within(o, hessian <- NA)), NA)


## egf_(expand|condense)_par ######
o <- egf_cache("egf-2.rds")
par <- o$tmb_out$env$last.par.best
len <- c(table(factor(names(par), levels = c("beta", "theta", "b"))))
epar <- egf_expand_par(o$tmb_out, par = par)
all.equal(epar, structure(c(par[1:4], theta = 0, par[-(1:4)]),
                             lengths = len + c(0L, 1L, 0L)))
cpar <- egf_condense_par(o$tmb_out, par = epar)
all.equal(cpar, structure(par, lengths = len))


## egf_get_sdreport ######
o <- egf_cache("egf-1.rds")
sdr <- o$sdreport
identical(egf_get_sdreport(o)$cov.fixed, sdr$cov.fixed)
o$sdreport <- NULL
assertWarning(all.equal(egf_get_sdreport(o)$cov.fixed, sdr$cov.fixed))


## egf_preprofile ######
## FIXME: Only testing a trivial case as no elements of 'beta' are mapped
o <- egf_cache("egf-1.rds")
l <- egf_preprofile(o, subset = seq_len(20L), top = c("log(r)", "log(c0)"))
is.list(l)
length(l) == 2L
identical(names(l), c("Y", "A"))
identical(unname(l$Y), Matrix(0, 20L, 2L))
identical(unname(l$A), sparseMatrix(i = seq_len(40L),
                                           j = rep.int(1:2, c(20L, 20L)),
                                           x = 1,
                                           dims = c(40L, 5L)))


## egf_cache ######
subdir <- tools::R_user_dir("epigrowthfit", "cache")
lf1 <- list.files(subdir)
file <- "test.rds"
a <- 1
x <- egf_cache(file, a)
y <- egf_cache(file, a <- 2)
lf2 <- list.files(subdir)
z <- egf_cache(file, clear = TRUE)
lf3 <- list.files(subdir)

identical(x, 1)
identical(y, 1)
identical(z, 0L)
identical(a, 1)
identical(setdiff(lf2, lf1), file)
identical(lf3, lf1)

