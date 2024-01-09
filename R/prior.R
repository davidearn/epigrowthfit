Normal <-
function(mu = 0, sigma = 1) {
	stopifnot(is.numeric(mu),
	          length(mu) > 0L,
	          is.finite(mu),
	          is.numeric(sigma),
	          length(sigma) > 0L,
	          is.finite(sigma),
	          sigma > 0)
	res <- list(family = "norm",
	            parameters = list(mu = as.double(mu), sigma = as.double(sigma)))
	class(res) <- "egf_prior"
	res
}

LKJ <-
function(eta = 1) {
	stopifnot(is.numeric(eta),
	          length(eta) > 0L,
	          is.finite(eta),
	          eta > 0)
	res <- list(family = "lkj",
	            parameters = list(eta = as.double(eta)))
	class(res) <- "egf_prior"
	res
}

Wishart <-
function(df, scale, tol = 1e-06) {
	stopifnot(is.numeric(tol),
	          length(tol) == 1L,
	          is.finite(tol),
	          tol >= 0)
	if (is.matrix(scale)) {
		scale <- list(scale)
	} else {
		stopifnot(is.list(scale))
	}
	for (i in seq_along(scale)) {
		stopifnot(is.numeric(scale[[i]]),
		          length(scale[[i]]) > 0L,
		          is.finite(scale[[i]]),
		          isSymmetric.matrix(scale[[i]]),
		          (e <- eigen(scale[[i]], symmetric = TRUE, only.values = TRUE)$values) > -tol * abs(e[1L]),
		          diag(scale[[i]]) > 0)
	}
	stopifnot(is.numeric(df),
	          length(df) > 0L,
	          is.finite(df),
	          rep.int(df, length(scale)) > rep.int(vapply(scale, nrow, 0L), length(df)) - 1L)
	scale <- lapply(scale, cov2theta)
	res <- list(family = "wishart",
	            parameters = list(df = as.double(df),
	                              scale = unname(scale)))
	class(res) <- "egf_prior"
	res
}

InverseWishart <-
function(df, scale, tol = 1e-06) {
	res <- Wishart(df = df, scale = scale, tol = tol)
	res$family <- "invwishart"
	res
}
