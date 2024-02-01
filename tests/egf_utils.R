attach(asNamespace("epigrowthfit"))
library(methods)
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## egf_sanitize_formula ################################################

l1 <- list(cbind(x, y) ~ 1,
           cbind(x, y) ~ g,
           cbind(x, y) ~ 1 + g,
           cbind(x, y) ~ (g),
           cbind(x, y) ~ g:h,
           cbind(x, y) ~ I(g + h),
           cbind(x, y) ~ I(g * h),
           cbind(x - 1, cumsum(y)) ~ g)
l2 <- list(~g,
           cbind(x, y) ~ g + h,
           cbind(x, y) ~ g * h,
           cbind(x, y) ~ 0 + g,
           cbind(x, y) ~ g - 1,
           cbind(x, y) ~ offset(h) + g,
           (cbind(x, y)) ~ g,
           cbind(x) ~ g,
           cbind(x, y, z) ~ g,
           rbind(x, y) ~ g) # i.e., anything other than 'cbind'

stopifnot(identical(lapply(l1, egf_sanitize_formula),
                    l1[c(1L, 2L, 2L, 2L, 5:8)]))
for (formula in l2)
	assertError(egf_sanitize_formula(formula))


## egf_sanitize_formula_parameters #####################################

model <- egf_model(curve = "exponential", family = "pois")
names_parameters <- egf_top(model, link = TRUE)

s <-
function(formula_parameters)
	egf_sanitize_formula_parameters(
		formula_parameters = formula_parameters,
		names_parameters = names_parameters,
		check_intercept = TRUE)

fp1 <- ~x * y + (z | g) + (zz | g/h)
l1 <- rep.int(list(simplify_terms(fp1)), 2L)
names(l1) <- c("log(r)", "log(c0)")

fp2 <- list(replace(fp1, 2:3, list(quote(log(r)), fp1[[2L]])))
l2 <- replace(l1, "log(c0)", list(~1))

fp3 <- c(fp2, list(log(c0) ~ x))
l3 <- replace(l2, "log(c0)", list(~x))

stopifnot(exprs = {
	identical(s(fp1), l1)
	identical(s(fp2), l2)
	identical(s(fp3), l3)
})
assertWarning(s(~0 + x))


## egf_make_frame ######################################################

model <- egf_model(curve = "exponential", family = "pois")

formula_ts <- cbind(day, count) ~ country
formula_windows <- cbind(left, right) ~ country
formula_parameters <- list(`log(r)`  = ~x1 + (1 | g1) + (1 | g1:g2),
                           `log(c0)` = ~(1 | g3))

data_ts <- data.frame(country = gl(6L, 11L),
                      day = seq.int(0, 10, by = 1),
                      count = rpois(11L, 100 * exp(0.04 * 0:10)))
data_windows <- data.frame(country = gl(3L, 2L),
                           left = rep.int(c(0, 5), 3L),
                           right = rep.int(c(5, 10), 3L),
                           x1 = c(5.00, 8.34, -0.57, -7.19, -9.71, 1.25),
                           x2 = rnorm(6L),
                           x3 = rnorm(6L),
                           g1 = c("a", "b", "b", "b", "b", "a"),
                           g2 = c("c", "d", "d", "d", "c", "c"),
                           g3 = c("f", "f", "e", "e", "e", "f"))

subset_ts <- quote(day > 0)
subset_windows <- quote(x1 < 0)
select_windows <- quote(.)

na_action_ts <- "pass"
na_action_windows <- "omit"

frame <- egf_make_frame(model = model,
                        formula_ts = formula_ts,
                        formula_windows = formula_windows,
                        formula_parameters = formula_parameters,
                        data_ts = data_ts,
                        data_windows = data_windows,
                        subset_ts = subset_ts,
                        subset_windows = subset_windows,
                        select_windows = select_windows,
                        na_action_ts = na_action_ts,
                        na_action_windows = na_action_windows)

stopifnot(exprs = {
	is.list(frame)
	length(frame) == 4L
	identical(names(frame), c("ts", "windows", "parameters", "extra"))
})

l1 <- frame[["ts"]]
l1.e <- data.frame(ts = gl(2L, 10L, labels = 2:3),
                   window = factor(rep.int(c(NA, 1, 2, NA, 3, NA),
                                           c(1L, 4L, 5L, 1L, 4L, 5L)),
                                   labels = sprintf("window_%d", 1:3)),
                   time = rep.int(seq.int(1, 10, by = 1), 2L),
                   x = data_ts[["count"]][c(NA, 14:22, NA, 25:33)])
attr(l1.e, "first") <- c(1L, 5L, 11L)
attr(l1.e, "last") <- c(5L, 10L, 15L)
stopifnot(identical(l1, l1.e))

l2 <- frame[["windows"]]
l2.e <- data.frame(ts = factor(c(2, 2, 3)),
                   window = gl(3L, 1L, labels = sprintf("window_%d", 1:3)),
                   start = c(1, 5, 1),
                   end = c(5, 10, 5))
stopifnot(identical(l2, l2.e))

l3 <- frame[["parameters"]]
stopifnot(exprs = {
	is.list(l3)
	length(l3) == 2L
	identical(names(l3), c("log(r)", "log(c0)"))
})

l31 <- l3[["log(r)"]]
l31.e <- droplevels(data_windows[3:5, c("x1", "g1", "g2"), drop = FALSE])
attr(l31.e, "terms") <- terms(formula_parameters[["log(r)"]])
row.names(l31.e) <- NULL
stopifnot(identical(l31, l31.e))

l32 <- l3[["log(c0)"]]
l32.e <- droplevels(data_windows[3:5, "g3", drop = FALSE])
attr(l32.e, "terms") <- terms(formula_parameters[["log(c0)"]])
row.names(l32.e) <- NULL
stopifnot(identical(l32, l32.e))

l4 <- frame[["extra"]]
l4.e <- droplevels(data_windows[3:5, c("country", "left", "right", "x2", "x3"), drop = FALSE])
row.names(l4.e) <- NULL
stopifnot(identical(l4, l4.e))


## egf_make_priors #####################################################

p1 <- Normal(mu = 0, sigma = 1)
p2 <- Normal(mu = 1, sigma = c(0.5, 1))
p3 <- Normal(mu = -1, sigma = 2)
p4 <- LKJ(eta = 1)

formula_priors <- list(foo(bar) ~ p1,
                       baz ~ p1,
                       beta ~ p1,
                       theta[[1L]] ~ p1,
                       theta[2:3] ~ p2,
                       theta[-(1:5)] ~ p3,
                       theta[replace(logical(6L), 4L, TRUE)] ~ p1,
                       Sigma ~ p4)

top <- list(names = c("foo(bar)", "baz"), family = "norm")

beta <- list(length = 4L, family = "norm")
theta <- list(length = 6L, family = "norm")
Sigma <- list(length = 1L, family = c("lkj", "wishart", "invwishart"),
              rows = 4L)

priors <- egf_make_priors(formula_priors = formula_priors,
                          top = top,
                          beta = beta,
                          theta = theta,
                          Sigma = Sigma)

p2.elt <-
function(i) {
	p2[["parameters"]][["sigma"]] <- p2[["parameters"]][["sigma"]][[i]]
	p2
}

stopifnot(exprs = {
	is.list(priors)
	length(priors) == 2L
	identical(names(priors), c("top", "bottom"))

	identical(priors[["top"]],
	          `names<-`(list(p1, p1), top[["names"]]))
	identical(priors[["bottom"]],
	          list(beta = list(p1, p1, p1, p1),
	               theta = list(p1, p2.elt(1L), p2.elt(2L), p1, NULL, p3),
	               Sigma = list(p4)))
})


## egf_make_X ##########################################################

formula <- ~x + y + z
data <- list(x = double(10L), y = gl(2L, 5L), z = rnorm(10L))
data <- model.frame(formula, data = data)

X1 <- egf_make_X(formula, data = data, sparse = FALSE)
X2 <- model.matrix(formula, data = data)
a <- attributes(X2)
X2 <- structure(X2[, -2L], assign = a[["assign"]][-2L], contrasts = a[["contrasts"]])
stopifnot(identical(X1, X2))

X1 <- egf_make_X(formula, data = data, sparse = TRUE)
X2 <- Matrix::sparse.model.matrix(formula, data = data)
a <- attributes(X2)
X2 <- structure(X2[, -2L], assign = a[["assign"]][-2L], contrasts = a[["contrasts"]])
stopifnot(identical(X1, X2))


## egf_make_Z ##########################################################

bar <- quote(x - 1 | f:g)
data <- list(x = rnorm(10L), f = gl(5L, 2L), g = gl(2L, 5L))
data <- model.frame(~x - 1 + f:g, data = data)
Z1 <- egf_make_Z(bar, data = data)
Z2 <- new("dgCMatrix",
          Dim = c(10L, 6L),
          Dimnames = list(as.character(1:10),
                          sprintf("(x | %s)", c("f1:g1", "f2:g1", "f3:g1", "f3:g2", "f4:g2", "f5:g2"))),
          p = c(0L, 2L, 4L, 5L, 6L, 8L, 10L),
          i = 0:9,
          x = data[["x"]])
attr(Z2, "assign") <- rep.int(1L, 6L)
attr(Z2, "level") <- gl(6L, 1L, labels = c("1:1", "2:1", "3:1", "3:2", "4:2", "5:2"))
stopifnot(identical(Z1, Z2))


## egf_combine_X #######################################################

fixed <- list(a = ~x, b = ~f)
data <- data.frame(x = 1:6, f = gl(2L, 3L))
l <- egf_combine_X(fixed = fixed,
                   X = lapply(fixed, egf_make_X, data = data, sparse = FALSE))

stopifnot(exprs = {
	is.list(l)
	length(l) == 3L
	identical(names(l), c("X", "effects", "contrasts"))
})

X <- cbind(1, 1:6, 1, rep(0:1, each = 3L))
dimnames(X) <- list(as.character(1:6),
                    c("(Intercept)", "x", "(Intercept)", "f2"))
effects <- data.frame(bottom = disambiguate(rep.int("beta", 4L)),
                      top = gl(2L, 2L, labels = c("a", "b")),
                      term = factor(c("(Intercept)", "x", "(Intercept)", "f")),
                      colname = colnames(X))
contrasts <- list(f = getOption("contrasts")[["unordered"]])

stopifnot(exprs = {
	identical(l[["X"]], X)
	identical(l[["effects"]], effects)
	identical(l[["contrasts"]], contrasts)
})


## egf_combine_Z #######################################################

random <- list(a = quote(x | f), b = quote(y | g))
data <- data.frame(x = 1:6, y = -1, f = gl(2L, 3L), g = gl(3L, 2L))
l <- egf_combine_Z(random = random,
                   Z = lapply(random, egf_make_Z, data = data))

stopifnot(exprs = {
	is.list(l)
	length(l) == 3L
	identical(names(l), c("Z", "effects", "contrasts"))
})

Z <- new("dgCMatrix",
         Dim = c(6L, 10L),
         Dimnames = list(as.character(1:6),
                         sprintf("(%s | %s)", rep.int(c("(Intercept)", "x", "y"), c(5L, 2L, 3L)), rep.int(c("f1", "f2", "g1", "g2", "g3"), 2L))),
         p = cumsum(c(0L, rep.int(c(3L, 2L, 3L, 2L), c(2L, 3L, 2L, 3L)))),
         i = rep.int(0:5, 4L),
         x = c(rep.int(1, 12L), 1:6, rep.int(-1, 6L)))

cov <- interaction(l[["effects"]][c("term", "group")],
                   drop = TRUE, lex.order = TRUE)
levels(cov) <- disambiguate(rep.int("Sigma", 4L))
vec <- interaction(l[["effects"]][c("term", "group", "level")],
                   drop = TRUE, lex.order = TRUE)
levels(vec) <- disambiguate(rep.int("u", 10L))

bottom <- disambiguate(rep.int("b", 10L))
top <- rep.int(factor(rep.int(c("a", "b"), c(2L, 3L))), 2L)
term <- factor(rep.int(c("(Intercept)", "x", "y"), c(5L, 2L, 3L)))
group <- rep.int(factor(rep.int(c("f", "g"), c(2L, 3L))), 2L)
level <- rep.int(factor(c(1:2, 1:3)), 2L)
colname <- colnames(Z)
effects <- data.frame(cov, vec, bottom, top, term, group, level, colname)

stopifnot(exprs = {
	identical(l[["Z"]], Z)
	identical(l[["effects"]], effects)
	is.null(l[["contrasts"]])
})
