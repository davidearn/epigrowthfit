devtools::load_all(".")
library("testthat")

test_that("get_names_top", {
  names_top_all <- get_names_top(NULL, link = FALSE)
  expect_type(names_top_all, "character")
  expect_gt(length(names_top_all), 0L)
  expect_false(anyNA(names_top_all))
  expect_named(names_top_all, NULL)

  model <- egf_model(
    curve = "richards",
    excess = TRUE,
    family = "nbinom",
    day_of_week = TRUE
  )
  names_top <- get_names_top(egf_model(), link = FALSE)
  expect_type(names_top, "character")
  expect_gt(length(names_top), 0L)
  expect_true(all(names_top %in% names_top_all))
})

test_that("string_*_link", {
  names_top_link0 <- get_names_top(NULL, link = FALSE)
  names_link <- string_get_link(names_top_link0)

  expect_equal(names_link, ifelse(names_top_link0 == "p", "logit", "log"))
  expect_equal(string_get_link("invalid string"), NA_character_)
  expect_equal(string_get_link("log(r)"), NA_character_)

  names_top_link1 <- string_add_link(names_top_link0)
  expect_equal(string_add_link(names_top_link0), sprintf("%s(%s)", names_link, names_top_link0))
  expect_equal(string_add_link("invalid string"), NA_character_)
  expect_equal(string_add_link("log(r)"), NA_character_)

  expect_equal(string_remove_link(names_top_link1), names_top_link0)
  expect_equal(string_remove_link("invalid string"), NA_character_)
  expect_equal(string_remove_link("r"), NA_character_)

  expect_equal(string_extract_link(names_top_link1), names_link)
  expect_equal(string_extract_link("invalid string"), NA_character_)
  expect_equal(string_extract_link("r"), NA_character_)

  expect_equal(names_top_link1, get_names_top(NULL, link = TRUE))
})

test_that("match_link", {
  expect_identical(match_link("identity"), identity)
  expect_identical(match_link("log"), log)
  expect_identical(match_link("logit"), function(p) qlogis(p), ignore_function_env = TRUE)
  expect_error(match_link("invalid string"))

  expect_identical(match_link("identity", inverse = TRUE), identity)
  expect_identical(match_link("log", inverse = TRUE), exp)
  expect_identical(match_link("logit", inverse = TRUE), function(q) plogis(q), ignore_function_env = TRUE)
  expect_error(match_link("invalid string", inverse = TRUE))
})

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
  expect_equal(row.names(o1), as.character(seq_len(20L)))
  expect_identical(attr(o1, "formula"), formula)
  expect_equal(attr(o1, "first"), c(1L, 5L, 11L))
  expect_equal(attr(o1, "last"), c(5L, 10L, 15L))
  expect_equal(o1$ts, gl(2L, 10L, labels = 2:3))
  expect_equal(o1$window, factor(rep.int(c(NA, 1, 2, NA, 3, NA), c(1L, 4L, 5L, 1L, 4L, 5L)), labels = sprintf("window_%d", 1:3)))
  expect_equal(o1$time, rep.int(1:10, 2L))
  expect_equal(o1$x, data$count[c(NA, 14:22, NA, 25:33)])

  o2 <- res$frame_windows
  expect_type(o2, "list")
  expect_s3_class(o2, "data.frame")
  expect_length(o2, 4L)
  expect_named(o2, c("ts", "window", "start", "end"))
  expect_equal(row.names(o2), as.character(seq_len(3L)))
  expect_identical(attr(o2, "formula"), formula_windows)
  expect_equal(o2$ts, factor(c(2, 2, 3)))
  expect_equal(o2$window, gl(3L, 1L, labels = sprintf("window_%d", 1:3)))
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
  expect_equal(row.names(o4), as.character(seq_len(3L)))
  expect_identical(attr(o4, "terms"), terms(formula_parameters$`log(r)`))
  expect_identical(o4, droplevels(data_windows[3:5, names(o4), drop = FALSE]), ignore_attr = c("row.names", "terms"))

  o5 <- res$frame_parameters$`log(c0)`
  expect_type(o5, "list")
  expect_s3_class(o5, "data.frame")
  expect_length(o5, 1L)
  expect_named(o5, "g3")
  expect_equal(row.names(o5), as.character(seq_len(3L)))
  expect_identical(attr(o5, "terms"), terms(formula_parameters$`log(c0)`))
  expect_identical(o5, droplevels(data_windows[3:5, names(o5), drop = FALSE]), ignore_attr = c("row.names", "terms"))

  o6 <- res$frame_append
  expect_type(o6, "list")
  expect_s3_class(o6, "data.frame")
  expect_length(o6, 5L)
  expect_named(o6, c("country", "left", "right", "x2", "x3"), ignore.order = TRUE)
  expect_equal(row.names(o6), as.character(seq_len(3L)))
  expect_identical(o6, droplevels(data_windows[3:5, names(o6), drop = FALSE]), ignore_attr = c("row.names", "terms"))
})

devtools::load_all(".")
test_that("egf_make_priors_top", {
  model <- egf_model()
  names_top <- get_names_top(model, link = TRUE)
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
  for (s in setdiff(names(priors), "log(r)")) {
    eval(bquote(expect_null(priors[[.(s)]])))
  }
})

test_that("egf_make_priors_bottom", {
  beta <- numeric(4L)
  theta <- numeric(10L)
  prior <- Normal(mu = 0, sigma = 1)
  formula_priors_bottom <- list(
    beta ~ prior,
    theta[[1L]] ~ prior,
    theta[2:3] ~ prior,
    theta[-c(1:3, 8:10)] ~ prior,
    theta[rep.int(c(FALSE, TRUE, FALSE), c(7L, 2L, 1L))] ~ prior
  )
  priors_bottom <- egf_make_priors_bottom(
    formula_priors_bottom = formula_priors_bottom,
    beta = beta,
    theta = theta
  )

  expect_type(priors_bottom, "list")
  expect_length(priors_bottom, sum(len))
  expect_named(priors_bottom, enum_dupl_string(rep.int(names(len), c(4L, 10L))))
  expect_identical(unname(priors_bottom[-sum(len)]), rep_len(list(prior), 13L))
  expect_null(priors_bottom[[14L]])
})

devtools::load_all(".")
library("testthat")

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
  expect_equal(attr(X1, "assign"), attr(mm1, "assign")[-2L])
  expect_identical(attr(X1, "contrasts"), attr(mm1, "contrasts"))

  mm2 <- sparse.model.matrix(formula, data = data)
  X2 <- egf_make_X(formula, data = data, sparse = TRUE)
  expect_identical(as.matrix(X2), as.matrix(mm2[, -2L]))
  expect_equal(X2@assign, mm2@assign[-2L])
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
  expect_equal(Z@Dim, c(10L, 6L))
  expect_identical(Z@Dimnames, list(NULL, sprintf("(x | %s)", c("f1:g1", "f2:g1", "f3:g1", "f3:g2", "f4:g2", "f5:g2"))))
  expect_equal(Z@x, data$x)
  expect_equal(Z@i, 0:9)
  expect_equal(Z@p, c(0L, cumsum(c(2L, 2L, 1L, 1L, 2L, 2L))))
  expect_equal(Z@assign, rep_len(1L, 6L))
  expect_error(Z@contrasts)
  expect_equal(Z@group, gl(6L, 1L, labels = c("1:1", "2:1", "3:1", "3:2", "4:2", "5:2")))
})

# test_that("egf_make_XZ_info", {
#
# })

# test_that("egf_make_tmb_data", {
#
# })
#
# test_that("egf_make_tmb_parameters", {
#
# })
#
# test_that("egf_make_tmb_args", {
#
# })
#
# test_that("egf_make_combined", {
#
# })
