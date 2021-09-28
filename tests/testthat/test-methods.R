example("egf", package = "epigrowthfit", local = TRUE, echo = FALSE)
object <- readRDS(system.file("exdata", "egf.rds", package = "epigrowthfit", mustWork = TRUE))

# test_that("confint", {
#
# })

test_that("fitted", {
  fo <- fitted(object, se = TRUE)
  p <- length(object$frame_parameters)
  n <- nrow(object$frame_windows)

  expect_type(fo, "list")
  expect_s3_class(fo, c("egf_fitted", "data.frame"), exact = TRUE)
  expect_length(fo, 5L)
  expect_named(fo, c("top", "ts", "window", "estimate", "se"))
  expect_identical(nrow(fo), p * n)
  expect_identical(fo$top, gl(p, n, labels = names(object$frame_parameters)))
  expect_identical(fo$ts, rep.int(object$frame_windows$ts, p))
  expect_identical(fo$window, rep.int(object$frame_windows$window, p))
  expect_type(fo$estimate, "double")
  expect_type(fo$se, "double")
  expect_identical(coef(object, se = TRUE), fo)

  level <- 0.95
  q <- qchisq(level, df = 1)
  cifo <- confint(fo, level = level)
  expect_type(cifo, "list")
  expect_s3_class(cifo, "data.frame")
  expect_length(cifo, 6L)
  expect_named(cifo, c(names(fo)[1:4], "lower", "upper"))
  expect_identical(nrow(cifo), nrow(fo))
  expect_identical(unclass(cifo[1:4]), unclass(fo[1:4]))
  expect_equal(cifo$lower, fo$estimate - sqrt(q) * fo$se)
  expect_equal(cifo$upper, fo$estimate + sqrt(q) * fo$se)
  expect_identical(attr(cifo, "level"), level)
})

test_that("fixef", {
  fo <- fixef(object)
  p <- length(object$frame_parameters)

  expect_type(fo, "list")
  expect_s3_class(fo, c("egf_fixef", "data.frame"), exact = TRUE)
  expect_identical(nrow(fo), p)
  expect_length(fo, 5L)
  expect_named(fo, c("bottom", "top", "term", "colname", "estimate"))
  expect_identical(fo$bottom, names(object$best)[grep("^beta\\[", names(object$best))])
  expect_identical(fo$top, gl(p, 1L, labels = names(object$frame_parameters)))
  expect_identical(fo$term, gl(1L, p, labels = "(Intercept)"))
  expect_identical(fo$colname, rep_len("(Intercept)", p))
  expect_identical(fo$estimate, unname(object$best)[grep("^beta\\[", names(object$best))])
})

test_that("getCall", {
  o1 <- list(call = call("egf.method"))
  class(o1) <- "egf"
  expect_identical(getCall(o1), call("egf"))
  class(o1) <- "egf_no_fit"
  expect_identical(getCall(o1), call("egf"))

  o2 <- list(call = call("simulate.method"))
  class(o2) <- "egf_model_simulate"
  expect_identical(getCall(o2), call("simulate"))
})

test_that("print", {
  capture.output({
    expect_condition(print(object), regexp = NA)
    expect_identical(print(object), object)
    expect_invisible(print(object))
  })
})

# test_that("predict", {
#
# })

# test_that("profile", {
#
# })

test_that("ranef", {
  ro <- ranef(object)
  p <- length(object$frame_parameters)
  n <- nrow(object$frame_windows)

  country <- object$frame_parameters[[1L]]$country
  country <- gl(nlevels(country), p, p * n, labels = levels(country))
  wave <- object$frame_parameters[[1L]]$wave
  wave <- gl(nlevels(wave), p * nlevels(country), p * n, labels = levels(wave))

  theta <- unname(object$best)[grep("^theta\\[", names(object$best))]
  R <- diag(p)
  R[upper.tri(R)] <- theta[-seq_len(p)]
  RTR <- t(R) %*% R
  diag_D <- exp(theta[seq_len(p)]) / sqrt(diag(RTR))
  RTR[] <- diag_D * RTR * rep(diag_D, each = p)
  dimnames(RTR) <- rep_len(list(names(object$frame_parameters)), 2L)

  expect_type(ro, "list")
  expect_s3_class(ro, c("egf_ranef", "data.frame"), exact = TRUE)
  expect_identical(nrow(ro), p * n)
  expect_length(ro, 10L)
  expect_named(ro, c("cov", "vec", "bottom", "top", "term", "group", "level", "colname", "mode", "sd"))
  expect_identical(ro$cov, gl(1L, p * n, labels = enum_dupl_string("Sigma")))
  expect_identical(ro$vec, gl(n, p, labels = enum_dupl_string(rep_len("u", n))))
  expect_identical(ro$bottom, names(object$best)[grep("^b\\[", names(object$best))])
  expect_identical(ro$top, gl(p, 1L, p * n, labels = names(object$frame_parameters)))
  expect_identical(ro$term, gl(1L, p * n, labels = "(Intercept)"))
  expect_identical(ro$group, gl(1L, p * n, labels = "country:wave"))
  expect_identical(ro$level, interaction(country, wave, drop = TRUE, sep = ":", lex.order = FALSE))
  expect_identical(ro$colname, sprintf("((Intercept) | country%s:wave%s)", as.character(country), as.character(wave)))
  expect_identical(ro$mode, unname(object$best)[grep("^b\\[", names(object$best))])
  expect_equal(ro$sd, rep.int(exp(theta[seq_len(p)]), n))
  expect_equal(attr(ro, "Sigma"), `names<-`(list(RTR), enum_dupl_string("Sigma")))
})

# test_that("simulate", {
#
# })

test_that("summary", {
  so <- summary(object)
  fo <- fitted(object)

  expect_type(so, "list")
  expect_s3_class(so, "egf_summary")
  expect_length(so, 4L)
  expect_named(so, c("convergence", "value", "gradient", "fitted"), ignore.order = TRUE)
  expect_identical(so$convergence, object$optimizer_out$convergence)
  expect_identical(so$value, object$value)
  expect_identical(so$gradient, object$gradient)

  expect_type(so$fitted, "double")
  expect_identical(dim(so$fitted), c(6L, nlevels(fo$top)))
  expect_identical(dimnames(so$fitted), list(names(summary(0)), levels(fo$top)))
  for (s in colnames(so$fitted)) {
    eval(bquote(expect_equal(so$fitted[, .(s)], c(summary(fo$estimate[fo$top == .(s)])))))
  }

  capture.output({
    expect_condition(print(so), regexp = NA)
    expect_identical(print(so), so)
    expect_invisible(print(so))
  })
})

test_that("vcov", {
  S <- vcov(object, full = TRUE, cor = FALSE)
  S0 <- vcov(object, full = FALSE, cor = FALSE)
  k <- grep("^beta\\[", names(object$best)[object$nonrandom])

  expect_type(S, "double")
  expect_s3_class(S, c("egf_vcov", "matrix", "array"), exact = TRUE)
  expect_identical(dim(S), rep_len(length(object$nonrandom), 2L))
  expect_identical(dimnames(S), rep_len(list(names(object$best)[object$nonrandom]), 2L))
  expect_equal(as.numeric(S), as.numeric(object$sdreport$cov.fixed))
  expect_equal(vcov(object, full = TRUE, cor = TRUE), cov2cor(S))
  expect_equal(unclass(S0), S[k, k, drop = FALSE])
})
