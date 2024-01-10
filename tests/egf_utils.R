library(epigrowthfit)
options(warn = 2L, error = recover)


## egf_sanitize_formula ######
l1 <- list(cbind(x, y) ~ 1,
           cbind(x, y) ~ g,
           cbind(x, y) ~ 1 + g,
           cbind(x, y) ~ (g),
           cbind(x, y) ~ g:h,
           cbind(x, y) ~ I(g + h),
           cbind(x, y) ~ I(g * h),
           cbind(x - 1, cumsum(y)) ~ g)
e1 <- function(x, y) {
    eval(bquote(identical(egf_sanitize_formula(.(x)), .(y))))
}
Map(e1, l1, l1[c(1L, 2L, 2L, 2L, 5:8)])

l2 <- list(~g,
           cbind(x, y) ~ g + h,
           cbind(x, y) ~ g * h,
           cbind(x, y) ~ 0 + g,
           cbind(x, y) ~ g - 1,
           cbind(x, y) ~ offset(h) + g,
           (cbind(x, y)) ~ g,
           cbind(x) ~ g,
           cbind(x, y, z) ~ g,
           not_cbind(x, y) ~ g)
e2 <- function(x) {
    eval(bquote(assertError(egf_sanitize_formula(.(x)))))
}
lapply(l2, e2)


## egf_sanitize_formula_parameters ######
model <- egf_model(curve = "exponential", excess = FALSE,
                   family = "pois", day_of_week = FALSE)
names_parameters <- egf_get_names_top(model, link = TRUE)
sanitize <- function(x) {
    egf_sanitize_formula_parameters(x, names_parameters = names_parameters,
                                    check_intercept = TRUE)
}

fp1 <- ~x * y + (z | g) + (zz | g/h)
l1 <- rep.int(list(simplify_terms(fp1)), 2L)
names(l1) <- c("log(r)", "log(c0)")
identical(sanitize(fp1), l1)

fp2 <- list(replace(fp1, 2:3, list(quote(log(r)), fp1[[2L]])))
l2 <- replace(l1, "log(c0)", list(~1))
identical(sanitize(fp2), l2, ignore_formula_env = TRUE)

fp3 <- c(fp2, list(log(c0) ~ x))
l3 <- replace(l2, "log(c0)", list(~x))
identical(sanitize(fp3), l3, ignore_formula_env = TRUE)

assertWarning(sanitize(~0 + x))


## egf_make_frame ######
model <- egf_model(curve = "exponential", family = "pois")
formula_ts <- cbind(day, count) ~ country
formula_windows <- cbind(left, right) ~ country
formula_parameters <- list(`log(r)`  = ~x1 + (1 | g1) + (1 | g1:g2),
                           `log(c0)` = ~(1 | g3))
data_ts <- data.frame(country = gl(6L, 11L),
                      day = seq.int(0, 10, by = 1),
                      count = rpois(11L, 100 * exp(0.04 * 0:10)))
data_windows <- data.frame(country = gl(3L, 2L),
                           left = c(0, 5, 0, 5, 0, 5),
                           right = c(5, 10, 5, 10, 5, 10),
                           x1 = c(5.00, 8.34, -0.57, -7.19, -9.71, 1.25),
                           x2 = rnorm(6L),
                           x3 = rnorm(6L),
                           g1 = c("a", "b", "b", "b", "b", "a"),
                           g2 = c("c", "d", "d", "d", "c", "c"),
                           g3 = c("f", "f", "e", "e", "e", "f"))
subset_ts <- quote(day > 0)
subset_windows <-  quote(x1 < 0)
na_action_ts <- "pass"
na_action_windows <- "omit"
append <- quote(.)

res <- egf_make_frame(model = model,
                      formula_ts = formula_ts,
                      formula_windows = formula_windows,
                      formula_parameters = formula_parameters,
                      data_ts = data_ts,
                      data_windows = data_windows,
                      subset_ts = subset_ts,
                      subset_windows = subset_windows,
                      na_action_ts = na_action_ts,
                      na_action_windows = na_action_windows,
                      append = append)
is.list(res)
length(res) == 4L
identical(names(res), c("ts", "windows", "parameters", "append"))

o1 <- res$ts
o1_expected <-
    data.frame(ts = gl(2L, 10L, labels = 2:3),
               window = factor(rep.int(c(NA, 1, 2, NA, 3, NA),
                                       c(1L, 4L, 5L, 1L, 4L, 5L)),
                               labels = sprintf("window_%d", 1:3)),
               time = rep.int(seq.int(1, 10, by = 1), 2L),
               x = data_ts$count[c(NA, 14:22, NA, 25:33)])
attr(o1_expected, "first") <- c(1L, 5L, 11L)
attr(o1_expected, "last") <- c(5L, 10L, 15L)
identical(o1, o1_expected)

o2 <- res$windows
o2_expected <-
    data.frame(ts = factor(c(2, 2, 3)),
               window = gl(3L, 1L, labels = sprintf("window_%d", 1:3)),
               start = c(1, 5, 1),
               end = c(5, 10, 5))
identical(o2, o2_expected)

o3 <- res$parameters
is.list(o3)
length(o3) == 2L
identical(names(o3), c("log(r)", "log(c0)"))

o31 <- o3$`log(r)`
o31_expected <- droplevels(data_windows[3:5, c("x1", "g1", "g2"), drop = FALSE])
attr(o31_expected, "terms") <- terms(formula_parameters$`log(r)`)
row.names(o31_expected) <- NULL
identical(o31, o31_expected)

o32 <- o3$`log(c0)`
o32_expected <- droplevels(data_windows[3:5, "g3", drop = FALSE])
attr(o32_expected, "terms") <- terms(formula_parameters$`log(c0)`)
row.names(o32_expected) <- NULL
identical(o32, o32_expected)

o4 <- res$append
o4_expected <- droplevels(data_windows[3:5, c("country", "left", "right", "x2", "x3"), drop = FALSE])
row.names(o4_expected) <- NULL
identical(o4, o4_expected)


## egf_make_priors ######
top <- list(names = c("foo(bar)", "baz"), family = "norm")
beta <- list(length = 4L, family = "norm")
theta <- list(length = 6L, family = "norm")
Sigma <- list(length = 1L, family = c("lkj", "wishart", "invwishart"),
              rows = 4L)

p1 <- Normal(mu = 0, sigma = 1)
p2 <- Normal(mu = 1, sigma = c(0.5, 1))
p3 <- Normal(mu = -1, sigma = 2)
p4 <- LKJ(eta = 1)
f <- function(i) {p2$parameters$sigma <- p2$parameters$sigma[[i]]; p2}

formula_priors <- list(foo(bar) ~ p1,
                       baz ~ p1,
                       beta ~ p1,
                       theta[[1L]] ~ p1,
                       theta[2:3] ~ p2,
                       theta[-(1:5)] ~ p3,
                       theta[replace(logical(6L), 4L, TRUE)] ~ p1,
                       Sigma ~ p4)
priors <- egf_make_priors(formula_priors = formula_priors,
                          top = top,
                          beta = beta,
                          theta = theta,
                          Sigma = Sigma)

is.list(priors)
length(priors) == 2L
identical(names(priors), c("top", "bottom"))
identical(priors$top, `names<-`(list(p1, p1), top$names))
identical(priors$bottom,
                 list(beta = list(p1, p1, p1, p1),
                      theta = list(p1, f(1L), f(2L), p1, NULL, p3),
                      Sigma = list(p4)))


## egf_make_X ######
formula <- ~x + y + z
data <- list(x = double(10L), y = gl(2L, 5L), z = rnorm(10L))
data <- model.frame(formula, data = data)

X1 <- egf_make_X(formula, data = data, sparse = FALSE)
X2 <- model.matrix(formula, data = data)
a <- attributes(X2)
X2 <- structure(X2[, -2L], assign = a$assign[-2L], contrasts = a$contrasts)
identical(X1, X2)

X1 <- egf_make_X(formula, data = data, sparse = TRUE)
X2 <- Matrix::sparse.model.matrix(formula, data = data)
a <- attributes(X2)
X2 <- structure(X2[, -2L], assign = a$assign[-2L], contrasts = a$contrasts)
identical(X1, X2)


## egf_make_Z ######
bar <- quote(x - 1 | f:g)
data <- list(x = rnorm(10L), f = gl(5L, 2L), g = gl(2L, 5L))
data <- model.frame(~x - 1 + f:g, data = data)
Z1 <- egf_make_Z(bar, data = data)
Z2 <- Matrix::sparseMatrix(i = seq_len(10L),
                           j = rep.int(seq_len(6L),
                                       c(2L, 2L, 1L, 1L, 2L, 2L)),
                           x = data$x,
                           dims = c(10L, 6L),
                           dimnames = list(as.character(seq_len(10L)),
                                           sprintf("(x | %s)", c("f1:g1", "f2:g1", "f3:g1", "f3:g2", "f4:g2", "f5:g2"))))
attr(Z2, "assign") <- rep.int(1L, 6L)
attr(Z2, "level") <- gl(6L, 1L, labels = c("1:1", "2:1", "3:1", "3:2", "4:2", "5:2"))
identical(Z1, Z2)


## egf_combine_X ######
fixed <- list(a = ~x, b = ~f)
data <- data.frame(x = 1:6, f = gl(2L, 3L))
l <- egf_combine_X(fixed = fixed,
                   X = lapply(fixed, egf_make_X, data = data,
                              sparse = FALSE))

is.list(l)
length(l) == 3L
identical(names(l), c("X", "effects", "contrasts"))

X <- cbind(1, 1:6, 1, rep.int(c(0, 1), c(3L, 3L)))
rownames(X) <- as.character(seq_len(6L))
colnames(X) <- c("(Intercept)", "x", "(Intercept)", "f2")
effects <- data.frame(bottom = disambiguate(rep.int("beta", 4L)),
                      top = gl(2L, 2L, labels = c("a", "b")),
                      term = factor(c("(Intercept)", "x",
                                      "(Intercept)", "f")),
                      colname = colnames(X))
contrasts <- list(f = getOption("contrasts")[["unordered"]])

identical(l$X, X)
identical(l$effects, effects)
identical(l$contrasts, contrasts)


## egf_combine_Z ######
random <- list(a = quote(x | f), b = quote(y | g))
data <- data.frame(x = 1:6, y = -1, f = gl(2L, 3L), g = gl(3L, 2L))
l <- egf_combine_Z(random = random,
                   Z = lapply(random, egf_make_Z, data = data))

is.list(l)
length(l) == 3L
identical(names(l), c("Z", "effects", "contrasts"))

Z <- cbind(c(1, 1, 1, 0, 0, 0),
           c(0, 0, 0, 1, 1, 1),
           c(1, 1, 0, 0, 0, 0),
           c(0, 0, 1, 1, 0, 0),
           c(0, 0, 0, 0, 1, 1),
           c(1, 2, 3, 0, 0, 0),
           c(0, 0, 0, 4, 5, 6),
           c(-1, -1, 0, 0, 0, 0),
           c(0, 0, -1, -1, 0, 0),
           c(0, 0, 0, 0, -1, -1))
rownames(Z) <- as.character(seq_len(6L))
colnames(Z) <- sprintf("(%s | %s)", rep.int(c("(Intercept)", "x", "y"), c(5L, 2L, 3L)), rep.int(c("f1", "f2", "g1", "g2", "g3"), 2L))
Z <- as(Z, "dgCMatrix")

cov <- interaction(l$effects[c("term", "group")],
                   drop = TRUE, lex.order = TRUE)
levels(cov) <- disambiguate(rep.int("Sigma", 4L))
vec <- interaction(l$effects[c("term", "group", "level")],
                   drop = TRUE, lex.order = TRUE)
levels(vec) <- disambiguate(rep.int("u", 10L))
bottom <- disambiguate(rep.int("b", 10L))
top <- rep.int(factor(rep.int(c("a", "b"), c(2L, 3L))), 2L)
term <- factor(rep.int(c("(Intercept)", "x", "y"), c(5L, 2L, 3L)))
group <- rep.int(factor(rep.int(c("f", "g"), c(2L, 3L))), 2L)
level <- rep.int(factor(c(seq_len(2L), seq_len(3L))), 2L)
colname <- colnames(Z)
effects <- data.frame(cov, vec, bottom, top, term, group, level, colname)

identical(l$Z, Z)
identical(l$effects, effects)
is.null(l$contrasts)

