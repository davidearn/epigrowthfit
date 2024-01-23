attach(asNamespace("epigrowthfit"))
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## negate ##############################################################

x <- quote(x)
log_x <- call("log", x)

minus_x <- call("-", x)
identical(negate(x), minus_x)
identical(negate(minus_x), x)

minus_log_x <- call("-", log_x)
identical(negate(log_x), minus_log_x)
identical(negate(minus_log_x), log_x)

identical(negate(1), call("-", 1))
identical(negate(call("-", 1)), 1)
identical(negate(-1), call("-", -1))
identical(negate(call("-", -1)), -1)


## (un)?split_terms ####################################################

identical(split_terms(quote(+1)), list(quote(1)))
identical(split_terms(quote(-1)), list(quote(-1)))

x <- quote(1 + a * b - b)
l <- list(1, quote(a * b), quote(-b))
identical(split_terms(x), l)
identical(unsplit_terms(l), x)

x <- quote(w + (x | f/g) + (y + z | h))
l <- list(quote(w), quote(x | f/g), quote(y + z | h))
identical(split_terms(x), l)
identical(unsplit_terms(l), x)

is.null(unsplit_terms(list()))


## split_effects #######################################################

x <- y ~ x + (1 | f) + (a + b | g)
l <- list(fixed = y ~ x, random = list(quote(1 | f), quote(a + b | g)))
identical(split_effects(x), l)


## split_interaction ###################################################

x <- quote(a:b:I(f:g):log(h))
l <- list(quote(a), quote(b), quote(I(f:g)), quote(log(h)))
identical(split_interaction(x), l)
assertError(split_interaction(list()))


## gsub_bar_plus #######################################################

x1 <- ~x + (1 | f) + (a + b | g)
x2 <- call("+", call("+", quote(x), quote(1 + f)), quote(a + b + g))
x2 <- as.formula(call("~", x2))
identical(gsub_bar_plus(x1), x2)


## simplify_terms ######################################################

x1 <- ~0 + x * y - y
y1 <- formula(terms(x1, simplify = TRUE))
identical(simplify_terms(x1), y1)

x2 <- ~0 + x * y - y + (1 | f/g) + (a | f) + (0 + b | f:g)
y2 <- y1
y2[[2L]] <- call("+",
                 call("+", y2[[2L]], quote((a | f))),
                 quote((b - 1 | f:g)))
identical(simplify_terms(x2), y2)
