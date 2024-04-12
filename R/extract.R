fixef.egf <-
function(object, ...) {
	par <- coef(object, random = TRUE, full = TRUE)
	nms <- c("bottom", "top", "term", "colname")
	data.frame(object$effects$beta[nms],
	           estimate = par[labels(par) == "beta"],
	           row.names = NULL,
	           stringsAsFactors = FALSE)
}

fixef.egf_no_fit <- fixef.egf

ranef.egf <-
function(object, build_cov = FALSE, ...) {
	stopifnot(egf_has_random(object),
	          is_true_or_false(build_cov))

	par <- coef(object, random = TRUE, full = TRUE)
	nms <- c("cov", "vec", "bottom", "top", "term", "group", "level", "colname")
	res <- data.frame(object$effects$b[nms],
	                  mode = par[labels(par) == "b"],
	                  row.names = NULL,
	                  stringsAsFactors = FALSE)
	if (!build_cov)
		return(res)

	d1 <- object$tmb_out$env$data$block_rows
	d2 <- object$tmb_out$env$data$block_cols
	p <- as.integer(choose(d1 + 1L, 2L))
	theta <- split(par[labels(par) == "theta"], rep.int(seq_along(p), p))
	Sigma <- lapply(theta, theta2cov)
	names(Sigma) <- levels(object$effects$b$cov)
	tt <- table(object$effects$b$top, object$effects$b$cov)
	for (i in seq_along(Sigma))
		dimnames(Sigma[[i]])[1:2] <- list(rownames(tt)[tt[, i] > 0L])
	attr(res, "Sigma") <- Sigma
	res
}

ranef.egf_no_fit <- ranef.egf

vcov.egf <-
function(object, ...)
	egf_get_sdreport(object)$cov.fixed

getCall.egf <-
function(x, ...) {
	call <- NextMethod("getCall")
	call[[1L]] <- quote(egf)
	call
}

getCall.egf_no_fit <- getCall.egf

model.frame.egf <-
function(formula,
         which = c("ts", "windows", "parameters", "extra", "combined"),
         full = FALSE,
         top = egf_top(formula),
         ...) {
	which <- match.arg(which)
	if (which == "combined") {
		res <- do.call(cbind, unname(formula$frame$parameters))
		res <- cbind(res, formula$frame$extra)
		res[duplicated(names(res))] <- NULL
		return(res)
	}
	res <- formula$frame[[which]]
	switch(which,
	       "ts" = if (full) res else res[!is.na(res$window), , drop = FALSE],
	       "windows" = res,
	       "parameters" = res[[match.arg(top)]],
	       "extra" = res)
}

model.frame.egf_no_fit <- model.frame.egf

model.matrix.egf <-
function(object,
         which = c("fixed", "random"),
         top = NULL,
         random = NULL,
         ...) {
	which <- match.arg(which)

	if (is.null(top)) {
		## Return the combined fixed or random effects design matrix
		## from object internals
		name <- switch(which,
		               fixed = if (object$control$sparse_X) "Xs" else "Xd",
		               random = "Z")
		res <- object$tmb_out$env$data[[name]]

		## Append 'contrasts' but not 'assign', which only makes sense
		## for submatrices
		attr(res, "contrasts") <- object$contrasts[[substr(name, 1L, 1L)]]
		return(res)
	}

	top <- match.arg(top, egf_top(object))
	frame <- model.frame(object, which = "parameters", top = top)
	l <- split_effects(formula(terms(frame))) # list(fixed = <formula>, random = <list of calls to `|`>)

	if (which == "fixed") {
		## Return parameter-specific fixed effects design matrix
		res <- egf_make_X(fixed = l$fixed,
		                  data = frame,
		                  sparse = object$control$sparse_X)
		return(res)
	}

	any_random_effects <- length(l$random) > 0L

	if (is.null(random)) {
		if (!any_random_effects) {
			## Return empty sparse matrix with correct number of rows
			res <- object$tmb_out$env$data$Z[, integer(0L), drop = FALSE]
			return(res)
		}

		## Return parameter-specific combined random effects design matrix
		Z <- lapply(l$random, egf_make_Z, data = frame)
		res <- do.call(cbind, Z)

		## Append 'contrasts' but not 'assign', which only makes sense
		## for submatrices
		contrasts <- unlist1(lapply(Z, attr, "contrasts"))
		contrasts[duplicated(names(contrasts))] <- NULL
		attr(res, "contrasts") <- contrasts
		return(res)
	}

	if (!any_random_effects)
		stop(gettextf("expected %s = %s : mixed effects formula for parameter '%s' does not contain random effects terms",
		              "random", "NULL", top),
		     domain = NA)

	if (!any(l$random == random))
		stop(gettextf("expected %s = %s or in %s",
		              "random", "NULL", deparse(as.call(c(list(quote(expression)), l$random)))),
		     domain = NA)

	## Return term-specific random effects design matrix
	egf_make_Z(random = random, data = frame)
}

model.matrix.egf_no_fit <- model.matrix.egf

terms.egf <-
function(x, top = egf_top(x), ...) {
	top <- match.arg(top)
	frame <- model.frame(x, which = "parameters", top = top)
	terms(frame)
}

terms.egf_no_fit <- terms.egf

formula.egf <-
function(x,
         top = egf_top(x),
         split = FALSE, ...) {
	top <- match.arg(top)
	res <- formula(terms(x, top = top))
	if (split) {
		l <- split_effects(res)
		res <- l$fixed
		attr(res, "random") <- lapply(l$random, function(x) call("(", x))
	}
	res
}

formula.egf_no_fit <- formula.egf

nobs.egf <-
function(object, ...) {
	mf <- model.frame(object, which = "ts", full = FALSE)
	sum(!is.na(mf$x))
}

nobs.egf_no_fit <- nobs.egf

df.residual.egf <-
function(object, ...)
	as.double(nobs(object)) - sum(!object$random)

df.residual.egf_no_fit <- df.residual.egf

logLik.egf <-
function(object, ...) {
	res <- -object$value
	attr(res, "df") <- sum(!object$random)
	attr(res, "nobs") <- nobs(object)
	class(res) <- "logLik"
	res
}

extractAIC.egf <-
function(fit, scale, k = 2, ...) {
	ll <- logLik(fit)
	edf <- attr(ll, "df")
	c(edf, -2 * as.double(ll) + k * edf)
}
