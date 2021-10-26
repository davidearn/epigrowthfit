test_that("egf_sanitize_formula", {
  l1 <- list(
    cbind(x, y) ~ 1,
    cbind(x, y) ~ g,
    cbind(x, y) ~ 1 + g,
    cbind(x, y) ~ (g),
    cbind(x, y) ~ g:h,
    cbind(x, y) ~ I(g + h),
    cbind(x, y) ~ I(g * h),
    cbind(x - 1, cumsum(y)) ~ g
  )
  e1 <- function(x, y) {
    eval(bquote(expect_identical(egf_sanitize_formula(.(x)), .(y))))
  }
  Map(e1, l1, l1[c(1L, 2L, 2L, 2L, 5:8)])

  l2 <- list(
    ~g,
    cbind(x, y) ~ g + h,
    cbind(x, y) ~ g * h,
    cbind(x, y) ~ 0 + g,
    cbind(x, y) ~ g - 1,
    cbind(x, y) ~ offset(h) + g,
    (cbind(x, y)) ~ g,
    cbind(x) ~ g,
    cbind(x, y, z) ~ g,
    not_cbind(x, y) ~ g
  )
  e2 <- function(x) {
    eval(bquote(expect_error(egf_sanitize_formula(.(x)))))
  }
  lapply(l2, e2)
})

test_that("egf_sanitize_formula_parameters", {
  model <- egf_model(curve = "exponential", family = "pois")
  names_parameters <- egf_get_names_top(model, link = TRUE)
  sanitize <- function(x) {
    egf_sanitize_formula_parameters(x, names_parameters = names_parameters, check_intercept = TRUE)
  }

  fp1 <- ~x * y + (z | g) + (zz | g/h)
  l1 <- rep.int(list(simplify_terms(fp1)), 2L)
  names(l1) <- c("log(r)", "log(c0)")
  expect_identical(sanitize(fp1), l1)

  fp2 <- list(replace(fp1, 2:3, list(quote(log(r)), fp1[[2L]])))
  l2 <- replace(l1, "log(c0)", list(~1))
  expect_identical(sanitize(fp2), l2, ignore_formula_env = TRUE)

  fp3 <- c(fp2, list(log(c0) ~ x))
  l3 <- replace(l2, "log(c0)", list(~x))
  expect_identical(sanitize(fp3), l3, ignore_formula_env = TRUE)

  expect_warning(sanitize(~0 + x))
})

test_that("egf_make_frame", {
  model <- egf_model(curve = "exponential", family = "pois")
  formula_ts <- cbind(day, count) ~ country
  formula_windows <- cbind(left, right) ~ country
  formula_parameters <- list(
    `log(r)`  = ~x1 + (1 | g1) + (1 | g1:g2),
    `log(c0)` = ~(1 | g3)
  )
  data_ts <- data.frame(
    country = gl(6L, 11L),
    day = seq.int(0, 10),
    count = rpois(11L, 100 * exp(0.04 * 0:10))
  )
  data_windows <- data.frame(
    country = gl(3L, 2L),
    left = c(0, 5, 0, 5, 0, 5),
    right = c(5, 10, 5, 10, 5, 10),
    x1 = c(5.00, 8.34, -0.57, -7.19, -9.71, 1.25),
    x2 = rnorm(6L),
    x3 = rnorm(6L),
    g1 = c("a", "b", "b", "b", "b", "a"),
    g2 = c("c", "d", "d", "d", "c", "c"),
    g3 = c("f", "f", "e", "e", "e", "f")
  )
  subset_ts <- quote(day > 0)
  subset_windows <-  quote(x1 < 0)
  na_action_ts <- "pass"
  na_action_windows <- "omit"
  append <- quote(.)

  res <- egf_make_frame(
    model = model,
    formula_ts = formula_ts,
    formula_windows = formula_windows,
    formula_parameters = formula_parameters,
    data_ts = data_ts,
    data_windows = data_windows,
    subset_ts = subset_ts,
    subset_windows = subset_windows,
    na_action_ts = na_action_ts,
    na_action_windows = na_action_windows,
    append = append
  )
  expect_type(res, "list")
  expect_length(res, 4L)
  expect_named(res, c("ts", "windows", "parameters", "append"), ignore.order = TRUE)

  o1 <- res$ts
  expect_type(o1, "list")
  expect_s3_class(o1, "data.frame")
  expect_length(o1, 4L)
  expect_named(o1, c("ts", "window", "time", "x"))
  expect_identical(row.names(o1), as.character(seq_len(20L)))
  expect_identical(attr(o1, "first"), c(1L, 5L, 11L))
  expect_identical(attr(o1, "last"), c(5L, 10L, 15L))
  expect_identical(o1$ts, gl(2L, 10L, labels = 2:3))
  expect_identical(o1$window, factor(rep.int(c(NA, 1, 2, NA, 3, NA), c(1L, 4L, 5L, 1L, 4L, 5L)), labels = sprintf("window_%d", 1:3)))
  expect_equal(o1$time, rep.int(1:10, 2L))
  expect_equal(o1$x, data_ts$count[c(NA, 14:22, NA, 25:33)])

  o2 <- res$windows
  expect_type(o2, "list")
  expect_s3_class(o2, "data.frame")
  expect_length(o2, 4L)
  expect_named(o2, c("ts", "window", "start", "end"))
  expect_identical(row.names(o2), as.character(seq_len(3L)))
  expect_identical(o2$ts, factor(c(2, 2, 3)))
  expect_identical(o2$window, gl(3L, 1L, labels = sprintf("window_%d", 1:3)))
  expect_equal(o2$start, c(1, 5, 1))
  expect_equal(o2$end, c(5, 10, 5))

  o3 <- res$parameters
  expect_type(o3, "list")
  expect_length(o3, 2L)
  expect_named(o3, c("log(r)", "log(c0)"))

  o4 <- o3$`log(r)`
  expect_type(o4, "list")
  expect_s3_class(o4, "data.frame")
  expect_length(o4, 3L)
  expect_named(o4, c("x1", "g1", "g2"), ignore.order = TRUE)
  expect_identical(row.names(o4), as.character(seq_len(3L)))
  expect_identical(attr(o4, "terms"), terms(formula_parameters$`log(r)`))
  expect_identical(o4, droplevels(data_windows[3:5, names(o4), drop = FALSE]), ignore_attr = c("row.names", "terms"))

  o5 <- o3$`log(c0)`
  expect_type(o5, "list")
  expect_s3_class(o5, "data.frame")
  expect_length(o5, 1L)
  expect_named(o5, "g3")
  expect_identical(row.names(o5), as.character(seq_len(3L)))
  expect_identical(attr(o5, "terms"), terms(formula_parameters$`log(c0)`))
  expect_identical(o5, droplevels(data_windows[3:5, names(o5), drop = FALSE]), ignore_attr = c("row.names", "terms"))

  o6 <- res$append
  expect_type(o6, "list")
  expect_s3_class(o6, "data.frame")
  expect_length(o6, 5L)
  expect_named(o6, c("country", "left", "right", "x2", "x3"), ignore.order = TRUE)
  expect_identical(row.names(o6), as.character(seq_len(3L)))
  expect_identical(o6, droplevels(data_windows[3:5, names(o6), drop = FALSE]), ignore_attr = c("row.names", "terms"))
})

test_that("egf_make_priors", {
  top <- list(names = c("foo(bar)", "baz"), family = "norm")
  beta <- list(length = 4L, family = "norm")
  theta <- list(length = 6L, family = "norm")
  Sigma <- list(length = 1L, family = c("lkj", "wishart", "invwishart"), rows = 4L)

  p1 <- Normal(mu = 0, sigma = 1)
  p2 <- Normal(mu = 1, sigma = c(0.5, 1))
  p3 <- Normal(mu = -1, sigma = 2)
  p4 <- LKJ(eta = 1)

  formula_priors <- list(
    foo(bar) ~ p1,
    baz ~ p1,
    beta ~ p1,
    theta[[1L]] ~ p1,
    theta[2:3] ~ p2,
    theta[-(1:5)] ~ p3,
    theta[replace(logical(6L), 4L, TRUE)] ~ p1,
    Sigma ~ p4
  )
  priors <- egf_make_priors(
    formula_priors = formula_priors,
    top = top,
    beta = beta,
    theta = theta,
    Sigma = Sigma
  )

  expect_type(priors, "list")
  expect_length(priors, 2L)
  expect_named(priors, c("top", "bottom"), ignore.order = TRUE)

  expect_type(priors$top, "list")
  expect_length(priors$top, 2L)
  expect_named(priors$top, top$names)
  expect_identical(unname(priors$top), list(p1, p1))

  expect_type(priors$bottom, "list")
  expect_length(priors$bottom, 3L)
  expect_named(priors$bottom, c("beta", "theta", "Sigma"))

  f <- function(i) {p2$parameters$sigma <- p2$parameters$sigma[[i]]; p2}
  expect_identical(priors$bottom$beta, list(p1, p1, p1, p1))
  expect_identical(priors$bottom$theta, list(p1, f(1L), f(2L), p1, NULL, p3))
  expect_identical(priors$bottom$Sigma, list(p4))
})

test_that("egf_make_X", {
  formula <- ~x + y + z
  data <- list(x = double(10L), y = gl(2L, 5L), z = rnorm(10L))
  data <- model.frame(formula, data = data)

  mm1 <- model.matrix(formula, data = data)
  X1 <- egf_make_X(formula, data = data, sparse = FALSE)
  expect_identical(X1, mm1[, -2L], ignore_attr = c("assign", "contrasts"))
  expect_identical(attr(X1, "assign"), attr(mm1, "assign")[-2L])
  expect_identical(attr(X1, "contrasts"), attr(mm1, "contrasts"))

  mm2 <- Matrix::sparse.model.matrix(formula, data = data)
  X2 <- egf_make_X(formula, data = data, sparse = TRUE)
  expect_identical(as.matrix(X2), as.matrix(mm2[, -2L]))
  expect_identical(X2@assign, mm2@assign[-2L])
  expect_identical(X2@contrasts, mm2@contrasts)
})

test_that("egf_make_Z", {
  bar <- quote(x - 1 | f:g)
  data <- list(
    x = rnorm(10L),
    f = gl(5L, 2L),
    g = gl(2L, 5L)
  )
  data <- model.frame(~x - 1 + f:g, data = data)
  Z <- egf_make_Z(bar, data = data)
  expect_identical(Z@Dim, c(10L, 6L))
  expect_identical(Z@Dimnames, list(NULL, sprintf("(x | %s)", c("f1:g1", "f2:g1", "f3:g1", "f3:g2", "f4:g2", "f5:g2"))))
  expect_identical(Z@x, data$x)
  expect_identical(Z@i, 0:9)
  expect_identical(Z@p, c(0L, cumsum(c(2L, 2L, 1L, 1L, 2L, 2L))))
  expect_identical(Z@assign, rep.int(1L, 6L))
  expect_error(Z@contrasts)
  expect_identical(Z@level, gl(6L, 1L, labels = c("1:1", "2:1", "3:1", "3:2", "4:2", "5:2")))
})

test_that("egf_combine_X", {
  fixed <- list(a = ~x, b = ~f)
  data <- data.frame(x = 1:6, f = gl(2L, 3L))
  X <- lapply(fixed, egf_make_X, data = data, sparse = FALSE)
  l <- egf_combine_X(fixed = fixed, X = X)

  expect_type(l, "list")
  expect_length(l, 3L)
  expect_named(l, c("X", "effects", "contrasts"), ignore.order = TRUE)

  expect_true(is.matrix(l$X))
  expect_identical(dim(l$X), c(6L, 4L))
  expect_identical(dimnames(l$X), list(as.character(1:6), c("(Intercept)", "x", "(Intercept)", "f2")))
  X0 <- cbind(1, 1:6, 1, rep.int(c(0, 1), c(3L, 3L)))
  expect_identical(unname(l$X), X0)

  expect_type(l$effects, "list")
  expect_s3_class(l$effects, "data.frame")
  expect_length(l$effects, 4L)
  expect_identical(row.names(l$effects), as.character(seq_len(4L)))
  expect_named(l$effects, c("bottom", "top", "term", "colname"))
  expect_identical(l$effects$bottom, disambiguate(rep.int("beta", 4L)))
  expect_identical(l$effects$top, gl(2L, 2L, labels = c("a", "b")))
  expect_identical(l$effects$term, factor(c("(Intercept)", "x", "(Intercept)", "f")))
  expect_identical(l$effects$colname, colnames(l$X))

  expect_equal(l$contrasts, list(f = getOption("contrasts")[["unordered"]]))
})

test_that("egf_combine_Z", {
  random <- list(a = quote(x | f), b = quote(y | g))
  data <- data.frame(x = 1:6, y = -1, f = gl(2L, 3L), g = gl(3L, 2L))
  Z <- lapply(random, egf_make_Z, data = data)
  l <- egf_combine_Z(random = random, Z = Z)

  expect_type(l, "list")
  expect_length(l, 3L)
  expect_named(l, c("Z", "effects", "contrasts"), ignore.order = TRUE)

  expect_type(l$Z, "S4")
  expect_s4_class(l$Z, "dgCMatrix")
  expect_identical(dim(l$Z), c(6L, 10L))
  expect_identical(dimnames(l$Z), list(NULL, sprintf("(%s | %s)", rep.int(c("(Intercept)", "x", "y"), c(5L, 2L, 3L)), rep.int(c("f1", "f2", "g1", "g2", "g3"), 2L))))
  Z0 <- cbind(
    c(1, 1, 1, 0, 0, 0),
    c(0, 0, 0, 1, 1, 1),
    c(1, 1, 0, 0, 0, 0),
    c(0, 0, 1, 1, 0, 0),
    c(0, 0, 0, 0, 1, 1),
    c(1, 2, 3, 0, 0, 0),
    c(0, 0, 0, 4, 5, 6),
    c(-1, -1, 0, 0, 0, 0),
    c(0, 0, -1, -1, 0, 0),
    c(0, 0, 0, 0, -1, -1)
  )
  expect_identical(unname(as(l$Z, "matrix")), Z0)

  expect_type(l$effects, "list")
  expect_s3_class(l$effects, "data.frame")
  expect_length(l$effects, 8L)
  expect_identical(row.names(l$effects), as.character(seq_len(10L)))
  expect_named(l$effects, c("cov", "vec", "bottom", "top", "term", "group", "level", "colname"))
  expect_identical(l$effects$cov, `levels<-`(interaction(l$effects[c("term", "group")], drop = TRUE, lex.order = TRUE), disambiguate(rep.int("Sigma", 4L))))
  expect_identical(l$effects$vec, `levels<-`(interaction(l$effects[c("term", "group", "level")], drop = TRUE, lex.order = TRUE), disambiguate(rep.int("u", 10L))))
  expect_identical(l$effects$bottom, disambiguate(rep.int("b", 10L)))
  expect_identical(l$effects$top, rep.int(factor(rep.int(c("a", "b"), c(2L, 3L))), 2L))
  expect_identical(l$effects$term, factor(rep.int(c("(Intercept)", "x", "y"), c(5L, 2L, 3L))))
  expect_identical(l$effects$group, rep.int(factor(rep.int(c("f", "g"), c(2L, 3L))), 2L))
  expect_identical(l$effects$level, rep.int(factor(c(seq_len(2L), seq_len(3L))), 2L))
  expect_identical(l$effects$colname, colnames(l$Z))
  expect_identical(do.call(order, unname(l$effects[c("cov", "vec", "top")])), seq_len(10L))

  expect_null(l$contrasts)
})
