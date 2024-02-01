library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)

.S3method("split", "formula",
          function(x, f, drop = FALSE, ...) {
          	l <- epigrowthfit:::split_effects(x)
          	r <- l[["fixed"]]
          	attr(r, "random") <- lapply(l[["random"]], function(x) call("(", x))
          	r
          })

o.1 <- egf_cache("egf-1.rds")
o.2 <- egf_cache("egf-2.rds")


## coef ################################################################

nms <- c("beta", "theta", "b")

o.2c <- coef(o.2, full = FALSE)
o.2c.e <- structure(as.double(o.2c),
                    names = rep.int(nms, c(2L, 2L, 40L)),
                    lengths = c(beta = 2L, theta = 2L, b = 40L),
                    map = list(beta = NULL, theta = c(1L, 2L, NA), b = NULL),
                    full = FALSE,
                    class = "egf_coef")
stopifnot(identical(o.2c, o.2c.e))

o.2cl <- as.list(o.2c)
o.2cl.e <- split(as.double(o.2c), factor(names(o.2c), levels = nms))
for (s in nms)
	attr(o.2cl.e[[s]], "map") <- attr(o.2c, "map")[[s]]
attr(o.2cl.e, "full") <- attr(o.2cl, "full")
stopifnot(identical(o.2cl, o.2cl.e))


## fixef ###############################################################

o.2c <- coef(o.2, full = TRUE)
o.2f <- fixef(o.2)
o.2f.e <- data.frame(bottom = epigrowthfit:::disambiguate(rep.int("beta", 2L)),
                     top = gl(2L, 1L, labels = c("log(r)", "log(c0)")),
                     term = gl(1L, 2L, labels = "(Intercept)"),
                     colname = rep.int("(Intercept)", 2L),
                     estimate = o.2c[1:2])
stopifnot(identical(o.2f, o.2f.e))


## ranef ###############################################################

o.2c <- coef(o.2, full = TRUE)
o.2r <- ranef(o.2, build_cov = TRUE)

o.2r.e <- data.frame(cov = gl(1L, 40L, labels = epigrowthfit:::disambiguate("Sigma")),
                     vec = gl(20L, 2L, labels = epigrowthfit:::disambiguate(rep.int("u", 20L))),
                     bottom = epigrowthfit:::disambiguate(rep.int("b", 40L)),
                     top = gl(2L, 1L, 40L, labels = c("log(r)", "log(c0)")),
                     term = gl(1L, 40L, labels = "(Intercept)"),
                     group = gl(1L, 40L, labels = "country:wave"),
                     level = gl(20L, 2L, labels = paste0(gl(10L, 1L, 20L, labels = LETTERS[1:10]), ":", gl(2L, 10L))),
                     colname = sprintf("((Intercept) | country%s:wave%s)", gl(10L, 2L, 40L, labels = LETTERS[1:10]), gl(2L, 20L)),
                     mode = o.2c[-(1:5)])

Sigma <- list(theta2cov(o.2c[3:5]))
names(Sigma) <- levels(o.2r[["cov"]])
dimnames(Sigma[[1L]])[1:2] <- list(levels(o.2r[["top"]]))
attr(o.2r.e, "Sigma") <- Sigma

stopifnot(identical(o.2r, o.2r.e))


## vcov ################################################################

stopifnot(identical(vcov(o.1), o.1[["sdreport"]][["cov.fixed"]]))


## getCall #############################################################

o <- list(call = call("egf.method"))
class(o) <- "egf"
stopifnot(identical(getCall(o), call("egf")))
class(o) <- "egf_no_fit"
stopifnot(identical(getCall(o), call("egf")))


## model.frame #########################################################

mf0 <- model.frame(o.1, which = "ts", full = FALSE)
mf1 <- model.frame(o.1, which = "ts", full = TRUE)
mf2 <- model.frame(o.1, which = "windows")
mf3 <- model.frame(o.1, which = "parameters", top = "log(r)")
mf4 <- model.frame(o.1, which = "extra")
mf5 <- model.frame(o.1, which = "combined")

frame <- o.1[["frame"]]

d5 <- do.call(cbind, unname(frame[["parameters"]]))
d5 <- cbind(d5, frame[["extra"]])
d5 <- d5[!duplicated(names(d5))]

stopifnot(exprs = {
	identical(mf0, frame[["ts"]][!is.na(frame[["ts"]][["window"]]), , drop = FALSE])
	identical(mf1, frame[["ts"]])
	identical(mf2, frame[["windows"]])
	identical(mf3, frame[["parameters"]][["log(r)"]])
	identical(mf4, frame[["extra"]])
	identical(mf5, d5)
})


## model.matrix ########################################################

X <- model.matrix(o.1, which =  "fixed")
Z <- model.matrix(o.1, which = "random")

X1 <- model.matrix(o.1, which =  "fixed", top = "log(r)")
Z1 <- model.matrix(o.1, which = "random", top = "log(r)")
Z2 <- model.matrix(o.1, which = "random", top = "log(r)",
                   random = quote(1 | country:wave))

o.1f <- fixef(o.1)
o.1r <- ranef(o.1)

stopifnot(exprs = {
	identical(X,
	          structure(o.1[["tmb_out"]][["env"]][["data"]][["Xd"]],
	                    contrasts = o.1[["contrasts"]][["X"]]))
	identical(Z,
	          structure(o.1[["tmb_out"]][["env"]][["data"]][["Z"]],
	                    contrasts = o.1[["contrasts"]][["Z"]]))
	identical(as.double(X1),
	          as.double(X[, o.1f[["top"]] == "log(r)", drop = FALSE]))
	identical(as.double(Z1),
	          as.double(Z[, o.1r[["top"]] == "log(r)", drop = FALSE]))
	identical(Z2, Z1)
})


## terms ###############################################################

stopifnot(identical(terms(o.1, top = "log(r)"),
                    terms(model.frame(o.1, which = "parameters", top = "log(r)"))))


## formula #############################################################

o.1f.0 <- formula(o.1, top = "log(r)", split = FALSE)
o.1f.1 <- formula(o.1, top = "log(r)", split = TRUE)
stopifnot(exprs = {
	identical(o.1f.0, formula(terms(o.1, top = "log(r)")))
	identical(o.1f.1, split(o.1f.0))
})


## nobs ################################################################

x <- model.frame(o.1, which = "ts", full = FALSE)[["x"]]
stopifnot(identical(nobs(o.1), sum(!is.na(x))))


## df.residual #########################################################

stopifnot(identical(df.residual(o.1),
                    as.double(nobs(o.1)) - sum(!o.1[["random"]])))


## logLik ##############################################################

stopifnot(identical(logLik(o.1),
                    structure(-o.1[["value"]],
                              df = sum(!o.1[["random"]]),
                              nobs = nobs(o.1),
                              class = "logLik")))


## extractAIC ##########################################################

ll <- logLik(o.1)
df <- attr(ll, "df")
stopifnot(identical(extractAIC(o.1), c(df, -2 * as.double(ll) + 2 * df)))
