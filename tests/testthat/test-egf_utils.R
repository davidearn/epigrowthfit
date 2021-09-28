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
  sanitize <- function(x) egf_sanitize_formula_parameters(x, model = model, ignore_intercept = FALSE)

  fp1 <- ~x * y + (z | g) + (zz | g/h)
  l1 <- rep_len(list(simplify_terms(fp1)), 2L)
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

test_that("egf_make_frames", {
  model <- egf_model(curve = "exponential", family = "pois")
  formula <- cbind(day, count) ~ country
  formula_windows <- cbind(left, right) ~ country
  formula_parameters <- list(
    `log(r)`  = ~x1 + (1 | g1) + (1 | g1:g2),
    `log(c0)` = ~(1 | g3)
  )
  data <- data.frame(
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
  subset <- quote(day > 0)
  subset_windows <-  quote(x1 < 0)
  na_action <- "pass"
  na_action_windows <- "omit"
  append <- quote(.)

  res <- egf_make_frames(
    model = model,
    formula = formula,
    formula_windows = formula_windows,
    formula_parameters = formula_parameters,
    data = data,
    data_windows = data_windows,
    subset = subset,
    subset_windows = subset_windows,
    na_action = na_action,
    na_action_windows = na_action_windows,
    append = append
  )
  expect_type(res, "list")
  expect_length(res, 4L)
  expect_named(res, c("frame", "frame_windows", "frame_parameters", "frame_append"), ignore.order = TRUE)

  o1 <- res$frame
  expect_type(o1, "list")
  expect_s3_class(o1, "data.frame")
  expect_length(o1, 4L)
  expect_named(o1, c("ts", "window", "time", "x"))
  expect_identical(row.names(o1), as.character(seq_len(20L)))
  expect_identical(attr(o1, "formula"), formula)
  expect_identical(attr(o1, "first"), c(1L, 5L, 11L))
  expect_identical(attr(o1, "last"), c(5L, 10L, 15L))
  expect_identical(o1$ts, gl(2L, 10L, labels = 2:3))
  expect_identical(o1$window, factor(rep.int(c(NA, 1, 2, NA, 3, NA), c(1L, 4L, 5L, 1L, 4L, 5L)), labels = sprintf("window_%d", 1:3)))
  expect_equal(o1$time, rep.int(1:10, 2L))
  expect_equal(o1$x, data$count[c(NA, 14:22, NA, 25:33)])

  o2 <- res$frame_windows
  expect_type(o2, "list")
  expect_s3_class(o2, "data.frame")
  expect_length(o2, 4L)
  expect_named(o2, c("ts", "window", "start", "end"))
  expect_identical(row.names(o2), as.character(seq_len(3L)))
  expect_identical(attr(o2, "formula"), formula_windows)
  expect_identical(o2$ts, factor(c(2, 2, 3)))
  expect_identical(o2$window, gl(3L, 1L, labels = sprintf("window_%d", 1:3)))
  expect_equal(o2$start, c(1, 5, 1))
  expect_equal(o2$end, c(5, 10, 5))

  o3 <- res$frame_parameters
  expect_type(o3, "list")
  expect_length(o3, 2L)
  expect_named(o3, c("log(r)", "log(c0)"))

  o4 <- res$frame_parameters$`log(r)`
  expect_type(o4, "list")
  expect_s3_class(o4, "data.frame")
  expect_length(o4, 3L)
  expect_named(o4, c("x1", "g1", "g2"), ignore.order = TRUE)
  expect_identical(row.names(o4), as.character(seq_len(3L)))
  expect_identical(attr(o4, "terms"), terms(formula_parameters$`log(r)`))
  expect_identical(o4, droplevels(data_windows[3:5, names(o4), drop = FALSE]), ignore_attr = c("row.names", "terms"))

  o5 <- res$frame_parameters$`log(c0)`
  expect_type(o5, "list")
  expect_s3_class(o5, "data.frame")
  expect_length(o5, 1L)
  expect_named(o5, "g3")
  expect_identical(row.names(o5), as.character(seq_len(3L)))
  expect_identical(attr(o5, "terms"), terms(formula_parameters$`log(c0)`))
  expect_identical(o5, droplevels(data_windows[3:5, names(o5), drop = FALSE]), ignore_attr = c("row.names", "terms"))

  o6 <- res$frame_append
  expect_type(o6, "list")
  expect_s3_class(o6, "data.frame")
  expect_length(o6, 5L)
  expect_named(o6, c("country", "left", "right", "x2", "x3"), ignore.order = TRUE)
  expect_identical(row.names(o6), as.character(seq_len(3L)))
  expect_identical(o6, droplevels(data_windows[3:5, names(o6), drop = FALSE]), ignore_attr = c("row.names", "terms"))
})

test_that("egf_make_priors_top", {
  model <- egf_model()
  names_top <- egf_get_names_top(model, link = TRUE)
  prior <- Normal(mu = 0, sigma = 1)
  formula_priors_top <- list(
    log(r) ~ prior
  )
  priors_top <- egf_make_priors_top(
    formula_priors_top = formula_priors_top,
    model = model
  )

  expect_type(priors_top, "list")
  expect_length(priors_top, length(names_top))
  expect_named(priors_top, names_top)
  expect_identical(priors_top[["log(r)"]], prior)
  for (s in setdiff(names(priors_top), "log(r)")) {
    eval(bquote(expect_null(priors_top[[.(s)]])))
  }
})

test_that("egf_make_priors_bottom", {
  p1 <- Normal(mu = 0, sigma = 1)
  p2 <- LKJ(eta = 1)
  formula_priors_bottom <- list(
    beta ~ p1,
    theta[[1L]] ~ p1,
    theta[2:3] ~ p1,
    theta[-c(1:3, 8:10)] ~ p1,
    theta[rep.int(c(FALSE, TRUE, FALSE), c(7L, 2L, 1L))] ~ p1,
    Sigma ~ p2
  )
  priors_bottom <- egf_make_priors_bottom(
    formula_priors_bottom = formula_priors_bottom,
    beta_size = 4L,
    theta_size = 10L,
    block_rows = 4L
  )

  expect_type(priors_bottom, "list")
  expect_length(priors_bottom, 3L)
  expect_named(priors_bottom, c("beta", "theta", "Sigma"))

  expect_type(priors_bottom$beta, "list")
  expect_length(priors_bottom$beta, 4L)
  expect_named(priors_bottom$beta, enum_dupl_string(rep_len("beta", 4L)))
  expect_identical(unname(priors_bottom$beta), rep_len(list(p1), 4L))

  expect_type(priors_bottom$theta, "list")
  expect_length(priors_bottom$theta, 10L)
  expect_named(priors_bottom$theta, enum_dupl_string(rep_len("theta", 10L)))
  expect_identical(unname(priors_bottom$theta)[-10L], rep_len(list(p1), 9L))
  expect_null(priors_bottom$theta[[10L]])

  expect_type(priors_bottom$Sigma, "list")
  expect_length(priors_bottom$Sigma, 1L)
  expect_named(priors_bottom$Sigma, enum_dupl_string("Sigma"))
  expect_identical(unname(priors_bottom$Sigma), list(p2))
})

test_that("egf_make_X", {
  formula <- ~x + y + z
  data <- list(
    x = rep_len(0, 10L),
    y = gl(2L, 5L),
    z = rnorm(10L)
  )
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
  expect_identical(Z@assign, rep_len(1L, 6L))
  expect_error(Z@contrasts)
  expect_identical(Z@index, gl(6L, 1L, labels = c("1:1", "2:1", "3:1", "3:2", "4:2", "5:2")))
})

test_that("egf_combine_X", {
  fixed <- list(a = ~x, b = ~f)
  data <- data.frame(x = 1:6, f = gl(2L, 3L))
  l <- lapply(fixed, egf_make_X, data = data, sparse = FALSE)
  X <- egf_combine_X(fixed = fixed, X = l)
  expect_type(X, "double")
  expect_true(is.matrix(X))
  expect_identical(dim(X), c(6L, 4L))
  expect_identical(dimnames(X), list(as.character(1:6), c("(Intercept)", "x", "(Intercept)", "f2")))
  X0 <- cbind(1, 1:6, 1, rep.int(c(0, 1), c(3L, 3L)))
  expect_identical(unname(X), X0, ignore_attr = "info")

  info <- attr(X, "info")
  expect_type(info, "list")
  expect_s3_class(info, "data.frame")
  expect_length(info, 4L)
  expect_identical(row.names(info), as.character(seq_len(4L)))
  expect_named(info, c("bottom", "top", "term", "colname"))
  expect_identical(info$bottom, enum_dupl_string(rep_len("beta", 4L)))
  expect_identical(info$top, gl(2L, 2L, labels = c("a", "b")))
  expect_identical(info$term, factor(c("(Intercept)", "x", "(Intercept)", "f")))
  expect_identical(info$colname, colnames(X))
})

test_that("egf_combine_Z", {
  random <- list(a = quote(x | f), b = quote(y | g))
  data <- data.frame(x = 1:6, y = -1, f = gl(2L, 3L), g = gl(3L, 2L))
  l <- lapply(random, egf_make_Z, data = data)
  Z <- egf_combine_Z(random = random, Z = l)
  expect_type(Z, "S4")
  expect_s4_class(Z, "dgCMatrix")
  expect_identical(dim(Z), c(6L, 10L))
  expect_identical(dimnames(Z), list(NULL, sprintf("(%s | %s)", rep.int(c("(Intercept)", "x", "y"), c(5L, 2L, 3L)), rep.int(c("f1", "f2", "g1", "g2", "g3"), 2L))))
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
  expect_identical(unname(as.matrix(Z)), Z0)

  info <- attr(Z, "info")
  expect_type(info, "list")
  expect_s3_class(info, "data.frame")
  expect_length(info, 8L)
  expect_identical(row.names(info), as.character(seq_len(10L)))
  expect_named(info, c("cov", "vec", "bottom", "top", "term", "group", "level", "colname"))
  expect_identical(info$cov, `levels<-`(interaction(info[c("term", "group")], drop = TRUE, lex.order = TRUE), enum_dupl_string(rep_len("Sigma", 4L))))
  expect_identical(info$vec, `levels<-`(interaction(info[c("term", "group", "level")], drop = TRUE, lex.order = TRUE), enum_dupl_string(rep_len("u", 10L))))
  expect_identical(info$bottom, enum_dupl_string(rep_len("b", 10L)))
  expect_identical(info$top, rep.int(factor(rep.int(c("a", "b"), c(2L, 3L))), 2L))
  expect_identical(info$term, factor(rep.int(c("(Intercept)", "x", "y"), c(5L, 2L, 3L))))
  expect_identical(info$group, rep.int(factor(rep.int(c("f", "g"), c(2L, 3L))), 2L))
  expect_identical(info$level, rep.int(factor(c(seq_len(2L), seq_len(3L))), 2L))
  expect_identical(info$colname, colnames(Z))
  expect_identical(do.call(order, unname(info[c("cov", "vec", "top")])), seq_len(10L))
})
