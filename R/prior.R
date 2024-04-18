Normal <-
function(mu = 0, sigma = 1) {
	stopifnot(is.numeric(mu),
	          length(mu) > 0L,
	          all(is.finite(range(mu))),
	          is.numeric(sigma),
	          length(sigma) > 0L,
	          all(is.finite(range(sigma))),
	          min(sigma) > 0)
	ans <- list(family = "norm",
	            parameters = list(mu = as.double(mu), sigma = as.double(sigma)))
	class(ans) <- "egf_prior"
	ans
}

LKJ <-
function(eta = 1) {
	stopifnot(is.numeric(eta),
	          length(eta) > 0L,
	          all(is.finite(range(eta))),
	          min(eta) > 0)
	ans <- list(family = "lkj",
	            parameters = list(eta = as.double(eta)))
	class(ans) <- "egf_prior"
	ans
}

Wishart <-
function(df, scale, tol = 1e-06) {
	if (is.matrix(scale))
		scale <- list(scale)
	stopifnot(is.numeric(df),
	          length(df) > 0L,
	          all(is.finite(range(df))),
	          is.list(scale),
	          is.numeric(tol),
	          length(tol) == 1L,
	          is.finite(tol),
	          tol >= 0)
	for (i in seq_along(scale))
	stopifnot(is.numeric(scale[[i]]),
	          length(scale[[i]]) > 0L,
	          all(is.finite(range(scale[[i]]))),
	          isSymmetric.matrix(scale[[i]]),
	          min(e <- eigen(scale[[i]], symmetric = TRUE, only.values = TRUE)[["values"]]) > -tol * abs(e[1L]),
	          min(diag(scale[[i]], names = FALSE)) > 0)
	stopifnot(all(rep.int(df, length(scale)) > rep.int(vapply(scale, nrow, 0L), length(df)) - 1L))
	scale <- unname(lapply(scale, cov2theta))
	ans <- list(family = "wishart",
	            parameters = list(df = as.double(df), scale = scale))
	class(ans) <- "egf_prior"
	ans
}

InverseWishart <-
function(df, scale, tol = 1e-06) {
	ans <- Wishart(df = df, scale = scale, tol = tol)
	ans[["family"]] <- "invwishart"
	ans
}
