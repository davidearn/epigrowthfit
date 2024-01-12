library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)


## coef ######
o <- egf_cache("egf-2.rds")
co <- coef(o, full = FALSE)
is.double(co)
co_expected <- structure(as.double(co),
                         names = rep.int(c("beta", "theta", "b"),
                                         c(2L, 2L, 40L)),
                         lengths = c(beta = 2L, theta = 2L, b = 40L),
                         map = list(beta = NULL, theta = c(1L, 2L, NA),
                                    b = NULL),
                         full = FALSE,
                         class = "egf_coef")
identical(co, co_expected)

lco <- as.list(co)
lco_expected <- split(as.double(co),
                      factor(names(co), levels = c("beta", "theta", "b")))
for (s in names(lco_expected)) {
    attr(lco_expected[[s]], "map") <- attr(co, "map")[[s]]
}
attr(lco_expected, "full") <- FALSE
identical(lco, lco_expected)


## fixef ######
o <- egf_cache("egf-2.rds")
co <- coef(o, full = TRUE)
fo <- fixef(o)
fo_expected <- data.frame(bottom = disambiguate(rep.int("beta", 2L)),
                          top = gl(2L, 1L, labels = c("log(r)", "log(c0)")),
                          term = gl(1L, 2L, labels = "(Intercept)"),
                          colname = rep.int("(Intercept)", 2L),
                          estimate = co[seq_len(2L)])
identical(fo, fo_expected)


## ranef ######
o <- egf_cache("egf-2.rds")
co <- coef(o, full = TRUE)
ro <- ranef(o, build_cov = TRUE)
ro_expected <-
    data.frame(cov = gl(1L, 40L, labels = disambiguate("Sigma")),
               vec = gl(20L, 2L, labels = disambiguate(rep.int("u", 20L))),
               bottom = disambiguate(rep.int("b", 40L)),
               top = gl(2L, 1L, 40L, labels = c("log(r)", "log(c0)")),
               term = gl(1L, 40L, labels = "(Intercept)"),
               group = gl(1L, 40L, labels = "country:wave"),
               level = gl(20L, 2L, labels = paste0(gl(10L, 1L, 20L, labels = LETTERS[1:10]), ":", gl(2L, 10L))),
               colname = sprintf("((Intercept) | country%s:wave%s)", gl(10L, 2L, 40L, labels = LETTERS[1:10]), gl(2L, 20L)),
               mode = co[-seq_len(5L)])
identical(ro, ro_expected, ignore_attr = "Sigma")
Sigma <- attr(ro, "Sigma")
Sigma_expected <- list(theta2cov(co[3:5]))
names(Sigma_expected) <- levels(ro$cov)
dimnames(Sigma_expected[[1L]])[1:2] <- list(levels(ro$top))
identical(Sigma, Sigma_expected)


## vcov ######
o <- egf_cache("egf-1.rds")
identical(vcov(o), o$sdreport$cov.fixed)


## getCall ######
o <- list(call = call("egf.method"))
class(o) <- "egf"
identical(getCall(o), call("egf"))
class(o) <- "egf_no_fit"
identical(getCall(o), call("egf"))


## model.frame ######
o <- egf_cache("egf-1.rds")

mf0 <- model.frame(o, which = "ts", full = FALSE)
identical(mf0, o$frame$ts[!is.na(o$frame$ts$window), , drop = FALSE])

mf1 <- model.frame(o, which = "ts", full = TRUE)
identical(mf1, o$frame$ts)

mf2 <- model.frame(o, which = "windows")
identical(mf2, o$frame$windows)

mf3 <- model.frame(o, which = "parameters", top = "log(r)")
identical(mf3, o$frame$parameters[["log(r)"]])

mf4 <- model.frame(o, which = "append")
identical(mf4, o$frame$append)

mf5 <- model.frame(o, which = "combined")
dd <- do.call(cbind, unname(o$frame$parameters))
dd <- cbind(dd, o$frame$append)
dd[duplicated(names(dd))] <- NULL
identical(mf5, dd)


## model.matrix ######
o <- egf_cache("egf-1.rds")
fo <- fixef(o)
ro <- ranef(o)

X <- model.matrix(o, which = "fixed")
Z <- model.matrix(o, which = "random")
identical(X, structure(o$tmb_out$env$data$Xd,
                              contrasts = o$contrasts$X))
identical(Z, structure(o$tmb_out$env$data$Z,
                              contrasts = o$contrasts$Z))

X1 <- model.matrix(o, which = "fixed", top = "log(r)")
identical(as.double(X1),
                 as.double(X[, fo$top == "log(r)", drop = FALSE]))

Z1 <- model.matrix(o, which = "random", top = "log(r)")
identical(as.double(Z1),
                 as.double(Z[, ro$top == "log(r)", drop = FALSE]))

Z11 <- model.matrix(o, which = "random", top = "log(r)",
                    random = (1 | country:wave))
identical(Z1, Z11)


## terms ######
o <- egf_cache("egf-1.rds")
identical(terms(o, top = "log(r)"),
                 terms(model.frame(o, which = "parameters", top="log(r)")))


## formula ######
split.formula <- function(x) {
    l <- split_effects(x)
    res <- l$fixed
    attr(res, "random") <- lapply(l$random, function(x) call("(", x))
    res
}
o <- egf_cache("egf-1.rds")
fo0 <- formula(o, top = "log(r)", split = FALSE)
fo1 <- formula(o, top = "log(r)", split = TRUE)
identical(fo0, formula(terms(o, top = "log(r)")))
identical(fo1, split(fo0))


## nobs ######
o <- egf_cache("egf-1.rds")
mf <- model.frame(o, which = "ts", full = FALSE)
identical(nobs(o), sum(!is.na(mf$x)))


## df.residual ######
o <- egf_cache("egf-1.rds")
identical(df.residual(o), as.double(nobs(o)) - sum(!o$random))


## logLik ######
o <- egf_cache("egf-1.rds")
identical(logLik(o), structure(-o$value,
                                      df = sum(!o$random),
                                      nobs = nobs(o),
                                      class = "logLik"))


## extractAIC ######
o <- egf_cache("egf-1.rds")
ll <- logLik(o)
edf <- attr(ll, "df")
identical(extractAIC(o), c(edf, -2 * as.double(ll) + 2 * edf))
