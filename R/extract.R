coef.egf <-
function(object, full = FALSE, ...) {
	stopifnot(is_true_or_false(full))
	if (full) {
		res <- egf_expand_par(object$tmb_out, par = object$best)
	} else {
		res <- object$best
		attr(res, "lengths") <- lengths(object$tmb_out$env$parameters)
	}
	map <- lapply(object$tmb_out$env$parameters, attr, "map")
	if (!egf_has_random(object)) {
		map[names(map) != "beta"] <- list(NULL)
	}
	for (i in seq_along(map)) {
		if (!is.null(map[[i]])) {
			map[[i]] <- as.integer(map[[i]]) + 1L
			map[[i]][map[[i]] == 0L] <- NA
		}
	}
	attr(res, "map") <- map
	attr(res, "full") <- full
	class(res) <- "egf_coef"
	res
}

coef.egf_no_fit <-
function(object, full = FALSE, ...) {
	object$best <- object$init
	coef.egf(object, full = full, ...)
}

print.egf_coef <-
function(x, ...) {
	y <- x
	attributes(x)[c("full", "lengths", "map", "class")] <- NULL
	NextMethod("print")
	invisible(y)
}

as.list.egf_coef <-
function(x, ...) {
	names(x) <- NULL
	len <- attr(x, "lengths")
	f <- rep.int(gl(length(len), 1L, labels = names(len)), len)
	res <- split(x, f)
	map <- attr(x, "map")
	for (s in names(res)) {
		attr(res[[s]], "map") <- map[[s]]
	}
	attr(res, "full") <- attr(x, "full")
	res
}

as.data.frame.egf_coef <- as.data.frame.vector

fixef.egf <-
function(object, ...) {
	par <- coef(object, full = TRUE)
	nms <- c("bottom", "top", "term", "colname")
	data.frame(object$effects$beta[nms],
	           estimate = par[names(par) == "beta"],
	           row.names = NULL,
	           stringsAsFactors = FALSE)
}

fixef.egf_no_fit <- fixef.egf

ranef.egf <-
function(object, build_cov = FALSE, ...) {
	stopifnot(egf_has_random(object),
	          is_true_or_false(build_cov))

	par <- coef(object, full = TRUE)
	nms <- c("cov", "vec", "bottom", "top", "term", "group", "level", "colname")
	res <- data.frame(object$effects$b[nms],
	                  mode = par[names(par) == "b"],
	                  row.names = NULL,
	                  stringsAsFactors = FALSE)
	if (!build_cov) {
		return(res)
	}

	d1 <- object$tmb_out$env$data$block_rows
	d2 <- object$tmb_out$env$data$block_cols
	p <- as.integer(choose(d1 + 1L, 2L))
	theta <- split(par[names(par) == "theta"], rep.int(seq_along(p), p))
	Sigma <- lapply(theta, theta2cov)
	names(Sigma) <- levels(object$effects$b$cov)
	tt <- table(object$effects$b$top, object$effects$b$cov)
	for (i in seq_along(Sigma)) {
		dimnames(Sigma[[i]])[1:2] <- list(rownames(tt)[tt[, i] > 0L])
	}
	attr(res, "Sigma") <- Sigma
	res
}

ranef.egf_no_fit <- ranef.egf

vcov.egf <-
function(object, ...) {
	egf_get_sdreport(object)$cov.fixed
}

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
	if (which == "ts") {
		stopifnot(is_true_or_false(full))
		if (full) {
			return(res)
		} else {
			return(res[!is.na(res$window), , drop = FALSE])
		}
	}
	if (which == "parameters") {
		top <- match.arg(top)
		return(res[[top]])
	}
	res
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

	if (!any_random_effects) {
		stop1("Expected 'random = NULL': mixed effects formula ",
		      "for parameter ", sQuote(top), " does not contain ",
		      "random effects terms.")
	}

	if (!any(l$random == random)) {
		stop("Expected 'random = NULL' or 'random' matching one of:\n\n",
		     paste0("  ", l$random, collapse = "\n"))
	}

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
