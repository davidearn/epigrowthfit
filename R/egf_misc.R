egf_top <-
function(object, ...)
	UseMethod("egf_top", object)

egf_top.default <-
function(object, link = TRUE, ...) {
	stopifnot(is.null(object))
	names_top <- c("r", "alpha", "c0", "tinfl", "K",
	               "p", "a", "b", "disp", paste0("w", 1:6))
	if (!link) {
		return(names_top)
	}
	egf_link_add(names_top)
}

egf_top.egf_model <-
function(object, link = TRUE, ...) {
	names_top <- switch(object$curve,
	                    exponential = c("r", "c0"),
	                    subexponential = c("alpha", "c0", "p"),
	                    gompertz = c("alpha", "tinfl", "K"),
	                    logistic = c("r", "tinfl", "K"),
	                    richards = c("r", "tinfl", "K", "a"))
	if (object$excess) {
		names_top <- c(names_top, "b")
	}
	if (object$family == "nbinom") {
		names_top <- c(names_top, "disp")
	}
	if (object$day_of_week > 0L) {
		names_top <- c(names_top, paste0("w", 1:6))
	}
	if (!link) {
		return(names_top)
	}
	egf_link_add(names_top)
}

egf_top.egf <-
function(object, link = TRUE, ...)
	egf_top(object$model, link = link)

egf_top.egf_no_fit <- egf_top.egf

egf_has_random <-
function(object) {
	stopifnot(inherits(object, c("egf", "egf_no_fit")))
	ncol(object$tmb_out$env$data$Z) > 0L
}

egf_has_converged <-
function(object, tol = 1) {
	stopifnot(inherits(object, "egf"))
	object$optimizer_out$convergence == 0L &&
	    is.finite(object$value) &&
	    all(is.finite(object$gradient)) &&
	    max(abs(object$gradient)) < tol &&
	    object$hessian
}

##' Convert between condensed and full parameter vectors
##'
##' Condense the full bottom level parameter vector \code{c(beta, theta, b)}
##' to the representation used by \pkg{TMB} (as in, e.g., \code{last.par.best}),
##' which excludes mapped elements, or do the inverse operation.
##'
##' @param obj
##'   A list returned by \code{\link[TMB]{MakeADFun}}.
##' @param par
##'   A numeric vector.
##'
##' @return
##' \code{egf_expand_par} returns \code{c(beta, theta, b)}.
##' \code{egf_condense_par} returns \code{c(cbeta, ctheta, cb)},
##' where \code{cbeta} is the condensed representation of \code{beta},
##' and so on.
##' Attribute \code{lengths} preserves the length of each segment.

egf_expand_par <-
function(obj, par) {
	l <- obj$env$parList(par[obj$env$lfixed()], par)
	if (ncol(obj$env$data$Z) == 0L) {
		l[names(l) != "beta"] <- list(double(0L))
	}
	len <- lengths(l)
	res <- unlist1(l)
	names(res) <- rep.int(names(l), len)
	attr(res, "lengths") <- len
	res
}

egf_condense_par <-
function(obj, par) {
	parameters <- obj$env$parameters
	if (ncol(obj$env$data$Z) == 0L) {
		parameters[names(parameters) != "beta"] <- list(double(0L))
	}
	f <-
	function(x) {
		if (is.null(map <- attr(x, "map"))) {
			res <- seq_along(x)
			attr(res, "n") <- length(x)
		} else {
			res <- match(seq_len(attr(x, "nlevels")) - 1L, map)
			attr(res, "n") <- length(map)
		}
		res
	}
	index <- lapply(parameters, f)
	len <- vapply(index, attr, 0L, "n")
	l <- split(par, rep.int(gl(length(len), 1L, labels = names(len)), len))
	l <- Map(`[`, l, index)
	len <- lengths(l)
	res <- unlist1(l)
	names(res) <- rep.int(names(l), len)
	attr(res, "lengths") <- len
	res
}

##' Extract TMB-generated covariance information
##'
##' A utility for extracting or, if necessary, computing covariance information,
##' reused by various methods.
##'
##' @param object
##'   An \code{"\link{egf}"} object.
##'
##' @return
##' An \code{"\link{sdreport}"} object.

egf_get_sdreport <-
function(object) {
	stopifnot(inherits(object, "egf"))
	res <- object$sdreport
	if (is.null(res)) {
		if (egf_has_random(object)) {
			warning(gettextf("computing Hessian matrix for random effects model; retry with %s to avoid recomputation",
			                 "object = update(object, se = TRUE)"),
			        domain = NA)
		}
		call <- quote(sdreport(object$tmb_out,
		                       par.fixed = object$best[!object$random],
		                       getReportCovariance = FALSE))
		res <- tryCatch(eval(call),
		                error = function(e) `[[<-`(e, "call", call))
	}
	if (inherits(res, "error")) {
		stop(res)
	}
	res
}

##' Patch TMB-generated functions
##'
##' Define wrapper functions on top of \code{\link[TMB]{MakeADFun}}-generated
##' functions \code{fn} and \code{gr}, so that function and gradient evaluations
##' can retry inner optimization using fallback methods in the event that the
##' default method (usually \code{\link[TMB]{newton}}) fails.
##'
##' @param fn,gr
##'   Functions to be patched, assumed to be the so-named elements
##'   of a \code{\link[TMB]{MakeADFun}}-generated list object.
##' @param inner_optimizer
##'   A list of \code{"\link{egf_inner_optimizer}"} objects
##'   specifying inner optimization methods to be tried in turn.
##'
##' @return
##' A function.

egf_patch_fn <-
function(fn, inner_optimizer) {
	e <- environment(fn)
	if (!exists(".egf_env", where = e, mode = "environment", inherits = FALSE)) {
		e$.egf_env <- new.env(parent = emptyenv())
	}
	e$.egf_env$fn <- fn
	e$.egf_env$inner_optimizer <- inner_optimizer

	## For 'R CMD check'
	last.par <- random <- inner.method <- inner.control <- .egf_env <- NULL
	pfn <-
	function(x = last.par[-random], ...) {
		oim <- inner.method
		oic <- inner.control
		on.exit({
			inner.method <<- oim
			inner.control <<- oic
		})
		for (io in .egf_env$inner_optimizer) {
			inner.method <<- io$method
			inner.control <<- io$control
			v <- .egf_env$fn(x, ...)
			if (is.numeric(v) && length(v) == 1L && is.finite(v)) {
				return(v)
			}
		}
		NaN # no warning to avoid duplication of 'optim' and 'nlminb' warnings
	}
	environment(pfn) <- e
	pfn
}

egf_patch_gr <-
function(gr, inner_optimizer) {
	e <- environment(gr)
	if (!exists(".egf_env", where = e, mode = "environment", inherits = FALSE)) {
		e$.egf_env <- new.env(parent = emptyenv())
	}
	e$.egf_env$gr <- gr
	e$.egf_env$inner_optimizer <- inner_optimizer

	# For 'R CMD check'
	last.par <- random <- inner.method <- inner.control <- .egf_env <- NULL
	pgr <-
	function(x = last.par[-random], ...) {
		oim <- inner.method
		oic <- inner.control
		on.exit({
			inner.method <<- oim
			inner.control <<- oic
		})
		n <- length(x)
		for (io in .egf_env$inner_optimizer) {
			inner.method <<- io$method
			inner.control <<- io$control
			v <- .egf_env$gr(x, ...)
			if (is.numeric(v) && length(v) == n && all(is.finite(v))) {
				return(v)
			}
		}
		warning(gettextf("unable to evaluate '%s', returning %f",
		                 "gr(x)", NaN),
		        domain = NA)
		NaN # warning because length 1 result is unexpected
	}
	environment(pgr) <- e
	pgr
}

egf_preprofile <-
function(object, subset, top) {
	if (object$control$profile) {
		stop(gettextf("fixed effects coefficients profiled out of likelihood; retry with %s",
		              "object = update(object, control = egf_control(profile = FALSE))"),
		     domain = NA)
	}

	Y <- object$tmb_out$env$data$Y
	Y <- Y[subset, top, drop = FALSE]
	X <- model.matrix(object, "fixed")
	X <- X[subset, , drop = FALSE]

	c0 <- coef(object, full = FALSE)
	len <- attr(c0, "lengths")
	map <- attr(c0, "map")$beta
	if (is.null(map)) {
		argna <- logical(len[["beta"]])
		fmap <- gl(len[["beta"]], 1L)
	} else {
		argna <- is.na(map)
		fmap <- factor(map[!argna])
	}

	c1 <- coef(object, full = TRUE)
	beta <- c1[names(c1) == "beta"]

	ftop <- factor(fixef(object)$top, levels = top)
	J <- as(ftop, "sparseMatrix")
	A <- KhatriRao(J, X)

	B <- KhatriRao(J, t(beta))
	Y <- Y + tcrossprod(X[, argna, drop = FALSE], B[, argna, drop = FALSE])

	A <- tcrossprod(A[, !argna, drop = FALSE], as(fmap, "sparseMatrix"))
	A <- cbind(A, Matrix(0, nrow(A), len[["theta"]]))

	if (!all(rowSums(abs(A)) > 0)) {
		stop("found population fitted value not depending on any fixed effects coefficient")
	}
	list(Y = Y, A = A)
}

egf_cache <-
function(file, object, topic = NULL, clear = FALSE, clear_all = FALSE, ...) {
	subdir <- R_user_dir("epigrowthfit", "cache")
	if (clear_all) {
		return(unlink(subdir, recursive = TRUE))
	}
	path <- file.path(subdir, file)
	if (clear) {
		return(unlink(path))
	}
	if (file.exists(path)) {
		return(readRDS(path))
	}
	if (missing(object)) {
		if (is.null(topic)) {
			name <- sub("-\\d+\\.rds$", "", file)
			pattern <- paste0("^", gsub("-", "\\\\W", name))
			pattern <- paste0(pattern, "$")
			hs <- help.search(pattern,
			                  ignore.case = FALSE,
			                  package = "epigrowthfit",
			                  fields = "name",
			                  types = "help",
			                  verbose = FALSE)
			topic <- hs$matches$Topic
			if (length(topic) != 1L) {
				stop(gettextf("resource not created because '%s' matches %d topics",
				              "file", length(topic)),
				     domain = NA)
			}
		}
		example(topic, character.only = TRUE, package = "epigrowthfit",
		        local = TRUE, echo = FALSE)
		if (!file.exists(path)) {
			stop(gettextf("examples for topic \"%s\" were sourced but resource was not created",
			              topic),
			     domain = NA)
		}
		return(readRDS(path))
	}
	if (!dir.exists(subdir)) {
		dir.create(subdir)
	}
	saveRDS(object, file = path, ...)
	object
}
