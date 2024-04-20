egf_top <-
function(object, ...)
	UseMethod("egf_top", object)

egf_top.default <-
function(object, link = TRUE, ...) {
	stopifnot(is.null(object))
	top <- c("r", "alpha", "c0", "tinfl", "K",
	         "p", "a", "b", "disp", paste0("w", 1:6))
	if (link) egf_link_add(top) else top
}

egf_top.egf_model <-
function(object, link = TRUE, ...) {
	top <- switch(object[["curve"]],
	              exponential = c("r", "c0"),
	              subexponential = c("alpha", "c0", "p"),
	              gompertz = c("alpha", "tinfl", "K"),
	              logistic = c("r", "tinfl", "K"),
	              richards = c("r", "tinfl", "K", "a"))
	if (object[["excess"]])
		top <- c(top, "b")
	if (object[["family"]] == "nbinom")
		top <- c(top, "disp")
	if (object[["day_of_week"]] > 0L)
		top <- c(top, paste0("w", 1:6))
	if (link) egf_link_add(top) else top
}

egf_top.egf <-
function(object, link = TRUE, ...)
	egf_top(object[["model"]], link = link, ...)

egf_top.egf_no_fit <- egf_top.egf

egf_has_random <-
function(object, check = TRUE) {
	if (check) stopifnot(inherits(object, c("egf", "egf_no_fit")))
	ncol(object[["tmb_out"]][["env"]][["data"]][["Z"]]) > 0L
}

egf_has_converged <-
function(object, check = TRUE, tol = 1) {
	if (check) stopifnot(inherits(object, "egf"))
	object[["optimizer_out"]][["convergence"]] == 0L &&
		is.finite(object[["value"]]) &&
		all(is.finite(range(object[["gradient"]]))) &&
		max(abs(range(object[["gradient"]]))) < tol &&
		object[["hessian"]]
}

##' Convert Condensed and Full Parameter Vectors
##'
##' Condense the full bottom level parameter vector \code{c(beta, theta, b)}
##' to the representation used by \pkg{TMB} (as in, e.g., \code{last.par.best}),
##' which excludes mapped elements, or do the inverse operation.
##'
##' @param obj
##'   a list returned by \code{\link[TMB]{MakeADFun}}.
##' @param par
##'   a numeric vector.
##'
##' @return
##' \code{egf_par_expand} returns \code{c(beta, theta, b)}.
##' \code{egf_par_condense} returns \code{c(cbeta, ctheta, cb)},
##' where \code{cbeta} is the condensed representation of \code{beta},
##' and so on.
##' Attribute \code{len} preserves the length of each segment.

egf_par_expand <-
function(obj, par) {
	names(par) <- NULL
	l <- obj[["env"]][["parList"]](par[obj[["env"]][["lfixed"]]()], par)
	if (ncol(obj[["env"]][["data"]][["Z"]]) == 0L)
		l[names(l) != "beta"] <- list(double(0L))
	ans <- unlist1(l)
	attr(ans, "len") <- lengths(l)
	ans
}

egf_par_condense <-
function(obj, par) {
	names(par) <- NULL
	parameters <- obj[["env"]][["parameters"]]
	f <-
	function(x) {
		if (is.null(map <- attr(x, "map"))) {
			ans <- seq_along(x)
			attr(ans, "n") <- length(x)
		}
		else {
			ans <- match(seq_len(attr(x, "nlevels")) - 1L, map)
			attr(ans, "n") <- length(map)
		}
		ans
	}
	index <- lapply(parameters, f)
	len <- vapply(index, attr, 0L, "n")
	l <- split(par, rep.int(gl(length(len), 1L, labels = names(len)), len))
	l <- Map(`[`, l, index)
	ans <- unlist1(l)
	attr(ans, "len") <- lengths(l)
	ans
}

##' Extract \pkg{TMB}-Generated \dQuote{Reports}
##'
##' Utilities for retrieving (or, if necessary, computing) the result
##' of a call to \code{report} or \code{\link[TMB]{sdreport}}.
##'
##' @param object an \code{\link{egf}} object.
##'
##' @return A list.

egf_report <-
function(object, check = TRUE) {
	if (check) stopifnot(inherits(object, "egf"))
	ans <- object[["tmb_out"]][[".__egf__."]][["report"]]
	if (is.null(ans)) {
		call <- quote(object[["tmb_out"]][["report"]](object[["best"]]))
		ans <- object[["tmb_out"]][[".__egf__."]][["report"]] <-
			tryCatch(eval(call), error = function(e) `[[<-`(e, "call", call))
	}
	if (inherits(ans, "error"))
		stop(ans)
	ans
}

egf_adreport <-
function(object, check = TRUE) {
	if (check) stopifnot(inherits(object, "egf"))
	ans <- object[["tmb_out"]][[".__egf__."]][["adreport"]]
	if (is.null(ans)) {
		if (egf_has_random(object, check = check))
			message("computing a Hessian matrix ...")
		call <- quote(sdreport(object[["tmb_out"]],
		                       par.fixed = object[["best"]][!object[["random"]]],
		                       getReportCovariance = FALSE))
		ans <- object[["tmb_out"]][[".__egf__."]][["adreport"]] <-
			tryCatch(eval(call), error = function(e) `[[<-`(e, "call", call))
	}
	if (inherits(ans, "error"))
		stop(ans)
	ans
}

##' Patch \pkg{TMB}-Generated Functions
##'
##' Define wrapper functions on top of \code{\link[TMB]{MakeADFun}}-generated
##' functions \code{fn} and \code{gr}, so that function and gradient
##' evaluations can retry inner optimization using fallback methods
##' in the event that the default method (usually \code{\link[TMB]{newton}})
##' fails.
##'
##' @param fn,gr
##'   functions to be patched, assumed to be the so-named elements
##'   of a list returned by \code{\link[TMB]{MakeADFun}}.
##' @param inner_optimizer
##'   a list of \code{\link{egf_inner_optimizer}} objects
##'   specifying inner optimization methods to be tried in turn.
##'
##' @return
##' A function.

egf_patch_fn <-
function(fn, inner_optimizer) {
	e <- environment(fn)
	if (!exists(".__egf__.", where = e, mode = "environment", inherits = FALSE))
		e[[".__egf__."]] <- new.env(parent = emptyenv())
	e[[".__egf__."]][["fn"]] <- fn
	e[[".__egf__."]][["inner_optimizer"]] <- inner_optimizer

	## For 'R CMD check'
	last.par <- random <- inner.method <- inner.control <- .__egf__. <- NULL
	pfn <-
	function(x = last.par[-random], ...) {
		oim <- inner.method
		oic <- inner.control
		on.exit({
			inner.method <<- oim
			inner.control <<- oic
		})
		for (io in .__egf__.[["inner_optimizer"]]) {
			inner.method <<- io[["method"]]
			inner.control <<- io[["control"]]
			v <- .__egf__.[["fn"]](x, ...)
			if (is.numeric(v) && length(v) == 1L && is.finite(v))
				return(v)
		}
		NaN # no warning to avoid duplication of 'optim' and 'nlminb' warnings
	}
	environment(pfn) <- e
	pfn
}

egf_patch_gr <-
function(gr, inner_optimizer) {
	e <- environment(gr)
	if (!exists(".__egf__.", where = e, mode = "environment", inherits = FALSE))
		e[[".__egf__."]] <- new.env(parent = emptyenv())
	e[[".__egf__."]][["gr"]] <- gr
	e[[".__egf__."]][["inner_optimizer"]] <- inner_optimizer

	# For 'R CMD check'
	last.par <- random <- inner.method <- inner.control <- .__egf__. <- NULL
	pgr <-
	function(x = last.par[-random], ...) {
		oim <- inner.method
		oic <- inner.control
		on.exit({
			inner.method <<- oim
			inner.control <<- oic
		})
		n <- length(x)
		for (io in .__egf__.[["inner_optimizer"]]) {
			inner.method <<- io[["method"]]
			inner.control <<- io[["control"]]
			v <- .__egf__.[["gr"]](x, ...)
			if (is.numeric(v) && length(v) == n && all(is.finite(range(v))))
				return(v)
		}
		warning(gettextf("unable to evaluate %s, returning %f",
		                 "gr(x)", NaN),
		        domain = NA)
		NaN # warning because length 1 result is unexpected
	}
	environment(pgr) <- e
	pgr
}

egf_cache <-
function(file, object, topic = NULL, clear = FALSE, clearAll = FALSE, ...) {
	subdir <- R_user_dir("epigrowthfit", "cache")
	if (clearAll)
		return(unlink(subdir, recursive = TRUE))
	path <- file.path(subdir, file)
	if (clear)
		return(unlink(path))
	if (file.exists(path))
		return(readRDS(path))
	if (missing(object)) {
		if (is.null(topic)) {
			name <- sub("-\\d+\\.rds$", "", file)
			pattern <- paste0("^", gsub("-", "\\\\W", name), "$")
			hs <- help.search(pattern,
			                  ignore.case = FALSE,
			                  package = "epigrowthfit",
			                  fields = "name",
			                  types = "help",
			                  verbose = FALSE)
			topic <- hs[["matches"]][["Topic"]]
			if (length(topic) != 1L)
				stop(gettextf("resource not created because '%s' matches %d topics",
				              "file", length(topic)),
				     domain = NA)
		}
		example(topic, character.only = TRUE, package = "epigrowthfit",
		        local = TRUE, echo = FALSE)
		if (file.exists(path))
			return(readRDS(path))
		else stop(gettextf("examples for topic \"%s\" were sourced but resource was not created",
		                   topic),
		          domain = NA)
	}
	if (!dir.exists(subdir))
		dir.create(subdir)
	saveRDS(object, file = path, ...)
	object
}
