fixef.egf <-
function(object, ...) {
	par <- coef(object, random = TRUE, full = TRUE)
	nms <- c("bottom", "top", "term", "colname")
	data.frame(object[["effects"]][["beta"]][nms],
	           value = par[labels(par) == "beta"],
	           row.names = NULL,
	           stringsAsFactors = FALSE)
}

fixef.egf_no_fit <- fixef.egf

ranef.egf <-
function(object, makeSigma = FALSE, ...) {
	stopifnot(egf_has_random(object), isTrueFalse(makeSigma))

	par <- coef(object, random = TRUE, full = TRUE)
	nms <- c("cov", "vec", "bottom", "top", "term", "group", "level", "colname")
	ans <- data.frame(object[["effects"]][["b"]][nms],
	                  mode = par[labels(par) == "b"],
	                  row.names = NULL,
	                  stringsAsFactors = FALSE)
	if (!makeSigma)
		return(ans)

	n <- object[["tmb_out"]][["env"]][["data"]][["block_rows"]]
	p <- n + ((n - 1L) * n) %/% 2L

	tt <- table(object[["effects"]][["b"]][["top"]],
	            object[["effects"]][["b"]][["cov"]])

	theta <- split(par[labels(par) == "theta"], rep.int(seq_along(p), p))
	Sigma <- vector("list", length(theta))
	names(Sigma) <- levels(object[["effects"]][["b"]][["cov"]])
	for (j in seq_along(theta)) {
		tmp <- theta2cov(theta[[j]])
		dimnames(tmp) <- rep.int(list(rownames(tt)[tt[, j] > 0L]), 2L)
		Sigma[[j]] <- tmp
	}
	attr(ans, "Sigma") <- Sigma
	ans
}

ranef.egf_no_fit <- ranef.egf

vcov.egf <-
function(object, ...)
	egf_get_sdreport(object)[["cov.fixed"]]

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
	ans <-
		if (which != "combined")
			formula[["frame"]][[which]]
		else {
			tmp <- do.call(cbind,
			               unname(c(formula[["frame"]][["parameters"]],
			                        list(formula[["frame"]][["extra"]]))))
			tmp[duplicated(names(tmp))] <- NULL
			tmp
		}
	switch(which,
	       "ts" = if (full) ans else ans[!is.na(ans[["window"]]), ],
	       "windows" = ans,
	       "parameters" = ans[[match.arg(top)]],
	       "extra" = ans,
	       "combined" = ans)
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
		               fixed = if (object[["control"]][["sparse_X"]])
		                       	"Xs"
		                       else "Xd",
		               random = "Z")
		ans <- object[["tmb_out"]][["env"]][["data"]][[name]]

		## Append 'contrasts' but not 'assign', which only makes sense
		## for submatrices
		attr(ans, "contrasts") <- object[["contrasts"]][[substr(name, 1L, 1L)]]
		return(ans)
	}

	top <- match.arg(top)
	frame <- model.frame(object, which = "parameters", top = top)
	l <- split_effects(formula(terms(frame)))

	if (which == "fixed") {
		## Return parameter-specific fixed effects design matrix
		ans <- egf_make_X(fixed = l[["fixed"]],
		                  data = frame,
		                  sparse = object[["control"]][["sparse_X"]])
		return(ans)
	}

	if (is.null(random)) {
		if (length(l[["random"]]) == 0L)
			## Return empty sparse matrix with correct number of rows
			ans <- object[["tmb_out"]][["env"]][["data"]][["Z"]][, integer(0L), drop = FALSE]
		else {
			## Return parameter-specific combined random effects design matrix
			Z <- lapply(l[["random"]], egf_make_Z, data = frame)
			ans <- do.call(cbind, Z)

			## Append 'contrasts' but not 'assign', which only makes sense
			## for submatrices
			contrasts <- unlist1(lapply(Z, attr, "contrasts"))
			contrasts[duplicated(names(contrasts))] <- NULL
			attr(ans, "contrasts") <- contrasts
		}
		return(ans)
	}

	stopifnot(is.call(random), identical(random[[1L]], quote(`|`)))

	if (length(l[["random"]]) == 0L)
		stop(gettextf("expected %s = %s : mixed effects model formula for parameter '%s' does not contain random effects terms",
		              "random", "NULL", top),
		     domain = NA)
	if (!any(match(l[["random"]], list(random), 0L)))
		stop(gettextf("expected %s = %s or in %s",
		              "random", "NULL", deparse(as.expression(l[["random"]]))),
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
function(x, top = egf_top(x), split = FALSE, ...) {
	top <- match.arg(top)
	ans <- formula(terms(x, top = top))
	if (split) {
		l <- split_effects(ans)
		ans <- l[["fixed"]]
		attr(ans, "random") <- lapply(l[["random"]], function(x) call("(", x))
	}
	ans
}

formula.egf_no_fit <- formula.egf

nobs.egf <-
function(object, ...) {
	mf <- model.frame(object, which = "ts", full = FALSE)
	sum(!is.na(mf[["x"]]))
}

nobs.egf_no_fit <- nobs.egf

df.residual.egf <-
function(object, ...)
	as.double(nobs(object)) - sum(!object[["random"]])

df.residual.egf_no_fit <- df.residual.egf

logLik.egf <-
function(object, ...) {
	ans <- -object[["value"]]
	attr(ans, "df") <- sum(!object[["random"]])
	attr(ans, "nobs") <- nobs(object)
	class(ans) <- "logLik"
	ans
}

extractAIC.egf <-
function(fit, scale, k = 2, ...) {
	ll <- logLik(fit)
	edf <- attr(ll, "df")
	c(edf, -2 * as.double(ll) + k * edf)
}
