##' Utilities for data validation and processing
##'
##' Functions used inside of \code{\link{egf}} to validate and process
##' user-supplied formulae and data. Most edge cases should be caught
##' and handled here, so that functions called downstream can operate
##' without excessive conditional logic.
##'
##' @return
##' \describe{
##' \item{egf_sanitize_formula}{
##'   A formula.
##' }
##' \item{egf_sanitize_formula_parameters}{
##'   A list of formulae.
##' }
##' \item{egf_make_frame}{
##'   A (recursive) list of data frames.
##' }
##' \item{egf_make_priors}{
##'   A (recursive) list of \code{"\link{egf_prior"}} objects.
##' }
##' }

egf_sanitize_formula <-
function(formula) {
	s <- deparse1(substitute(formula))
	tt <- terms(formula, simplify = TRUE)
	a <- attributes(tt)
	if (!(a$response > 0L &&
	      is.call(lhs <- a$variables[[1L + a$response]]) &&
	      lhs[[1L]] == "cbind" &&
	      length(lhs) == 3L))
		stop(gettextf("left hand side of '%s' must be a call to '%s' with %d arguments",
		             s, "cbind", 2L),
		     domain = NA)
	if (!(a$intercept == 1L &&
	      is.null(a$offset) &&
	      length(a$term.labels) < 2L))
		stop(gettextf("right hand side of '%s' must be %d or have exactly %d term",
		              s, 1L, 1L),
		     domain = NA)
	formula(tt)
}

egf_sanitize_formula_parameters <-
function(formula_parameters, names_parameters, check_intercept) {
	recycled <- inherits(formula_parameters, "formula")
	if (recycled) {
		formula_parameters <- simplify_terms(formula_parameters)
		p <- length(names_parameters)
		formula_parameters <- rep.int(list(formula_parameters), p)
		names(formula_parameters) <- names_parameters
	}
	else {
		l <- lapply(formula_parameters, `[[`, 2L)
		names(formula_parameters) <- vapply(l, deparse, "")
		m <- match(names(formula_parameters), names_parameters, 0L)
		if (min(1L, m) == 0L)
			stop(gettextf("%s[[%d]][[2L]] is not an element of %s",
			              "formula_parameters", which.min(m), deparse(as.call(c(list(quote(expression)), l)))),
			     domain = NA)
		g <- function(x) simplify_terms(x[-2L])
		formula_parameters <- lapply(formula_parameters, g)
		need <- setdiff(names_parameters, names(formula_parameters))
		formula_parameters[need] <- list(`environment<-`(~1, globalenv()))
		formula_parameters <- formula_parameters[names_parameters]
	}
	if (check_intercept) {
		check <-
		function(x) {
			a <- attributes(terms(split_effects(x)$fixed))
			a$intercept == 1L || length(a$term.labels) == 0L
		}
		ok <-
			if (recycled)
				check(formula_parameters[[1L]])
			else all(vapply(formula_parameters, check, FALSE))
		if (!ok)
			warning(gettextf("default initial values of fixed effects coefficients not reliable for zero-intercept models; consider setting '%s' argument '%s'",
			                 "egf", "init"),
			        domain = NA)
	}
	formula_parameters
}

egf_make_frame <-
function(model,
         formula_ts, formula_windows, formula_parameters,
         data_ts, data_windows,
         subset_ts, subset_windows, select_windows,
         na_action_ts, na_action_windows) {
	## Reused for 'frame_ts' and 'frame_windows'
	make_frame <-
	function(formula, data, subset, na.action, drop.unused.levels) {
		group <- length(attr(terms(formula), "term.labels")) == 1L
		tt <- c(if (group) list(formula[[3L]]), as.list(formula[[2L]])[-1L])
		names_original <- c(if (!group) "", vapply(tt, deparse1, ""))
		tt <- Map(call, c(if (group) "as.factor", "identity", "identity"), tt,
		          USE.NAMES = FALSE)
		cl <- match.call()
		cl[[1L]] <- quote(model.frame)
		cl$formula <- as.formula(call("~", unsplit_terms(tt)),
		                         env = environment(formula))
		cl$subset <- subset
		mf <- eval(cl, parent.frame())
		if (!group)
			mf <- data.frame(gl(1L, nrow(mf)), mf)
		attr(mf, "names_original") <- names_original
		attr(mf, "terms") <- NULL
		mf
	}


	### Construction step ######################################################

	### Time series stuff

	frame_ts <- make_frame(formula = formula_ts,
	                       data = data_ts,
	                       subset = subset_ts,
	                       na.action = na.pass,
	                       drop.unused.levels = TRUE)
	names(frame_ts) <- c("ts", "time", "x")
	frame_ts <- frame_ts[!is.na(frame_ts$ts), , drop = FALSE]
	if (nrow(frame_ts) == 0L)
		stop(gettextf("data frame constructed from '%s' has zero rows",
		              "*_ts"),
		     domain = NA)

	### Fitting window stuff

	frame_windows <- make_frame(formula = formula_windows,
	                            data = data_windows,
	                            subset = subset_windows,
	                            na.action = switch(na_action_windows,
	                                               fail = na.fail,
	                                               omit = na.pass),
	                            drop.unused.levels = TRUE)
	names(frame_windows) <- c("ts", "start", "end")
	if ((N <- nrow(frame_windows)) == 0L)
		stop(gettextf("data frame constructed from '%s' has zero rows",
		              "*_windows"),
		     domain = NA)

	### Mixed effects stuff

	## Build model frames
	cl <- call("model.frame",
	           formula = NULL,
	           data = quote(data_windows),
	           subset = subset_windows,
	           na.action = switch(na_action_windows,
	                              fail = quote(na.fail),
	                              omit = quote(na.pass)),
	           drop.unused.levels = TRUE)
	frame_parameters <- vector("list", length(formula_parameters))
	names(frame_parameters) <- names(formula_parameters)
	for (i in seq_along(frame_parameters)) {
		cl$formula <- gsub_bar_plus(formula_parameters[[i]])
		frame_parameters[[i]] <- eval(cl)
	}

	## Replace, e.g., 'model.frame(~1, data)', which is an empty data frame
	## with 0 length and 0 rows regardless of 'data', with a data frame with
	## 0 length and the _correct_ number of rows
	len <- lengths(frame_parameters)
	frame_parameters[len == 0L] <- list(data.frame(row.names = seq_len(N)))

	## Test model frames for rowwise correspondence with 'frame_windows'
	if (any(vapply(frame_parameters, nrow, 0L) != N))
		stop(gettextf("data frames constructed from '%s' and '%s' do not have a common number of rows",
		              "formula_windows", "formula_parameters"),
		     domain = NA)

	### Appended stuff

	frame_extra <-
	if (!is.null(select_windows) && is.data.frame(data_windows)) {
		i <- egf_eval_subset(subset_windows, data_windows,
		                     environment(formula_windows))
		j <-
			if (select_windows == ".")
				setdiff(names(data_windows),
				        unlist1(lapply(frame_parameters, names)))
			else egf_eval_select(select_windows, data_windows, baseenv())
		data_windows[i, j]
	}
	else data.frame(row.names = seq_len(N))


	### Validation step ########################################################

	### Time series stuff

	nf1 <- as.list(attr(frame_ts, "names_original"))
	names(nf1) <- names(frame_ts)

	if (inherits(frame_ts$time, c("Date", "POSIXct", "POSIXlt")))
		frame_ts$time <- julian(frame_ts$time)
	else if (!is.numeric(frame_ts$time))
		stop(gettextf("'%s' is not a numeric, Date, or POSIXt vector", nf1$time),
		     domain = NA)
	if (!is.numeric(frame_ts$x))
		stop(gettextf("'%s' is not a numeric vector", nf1$x),
		     domain = NA)

	### Fitting window stuff

	nf2 <- as.list(attr(frame_windows, "names_original"))
	names(nf2) <- names(frame_windows)
	for (s in c("start", "end"))
		if (inherits(frame_windows[[s]], c("Date", "POSIXct", "POSIXlt")))
			frame_windows[[s]] <- julian(frame_windows[[s]])
		else if (!is.numeric(frame_windows[[s]]))
			stop(gettextf("'%s' is not a numeric, Date, or POSIXt vector", nf2[[s]]),
			     domain = NA)

	### Mixed effects stuff

	get_names_bar_lhs <-
	function(formula) {
		bars <- split_effects(formula)$random
		if (length(bars) == 0L)
			return(character(0L))
		f <-
		function(x) {
			vars <- attr(terms(as.formula(call("~", x[[2L]]))), "variables")
			vapply(vars, deparse1, "")[-1L]
		}
		unique(unlist1(lapply(bars, f)))
	}
	names_bar_lhs <- lapply(formula_parameters, get_names_bar_lhs)
	check_ok_bar_lhs <-
	function(data, names) {
		f <-
		function(x)
			is.double(x) || is.integer(x) || is.logical(x)
		all(vapply(data[names], f, FALSE))
	}
	if (!all(mapply(check_ok_bar_lhs, data = frame_parameters, names = names_bar_lhs)))
		stop(gettextf("variables on left hand side of '%s' in '%s' must be of double, integer, or logical type",
		              "|", "formula_parameters"),
		     domain = NA)


	### Subsetting step ########################################################

	## Discard fitting windows without observations on all
	## 'formula_windows' and 'formula_parameters' variables
	cc <- do.call(complete.cases,
	              c(list(frame_windows), frame_parameters[len > 0L]))
	if (!all(cc)) {
		frame_windows <- frame_windows[cc, , drop = FALSE]
		frame_parameters <- lapply(frame_parameters, `[`, cc, , drop = FALSE)
		frame_extra <- frame_extra[cc, , drop = FALSE]
		frame_windows$ts <- droplevels(frame_windows$ts)
	}

	## Discard time series without corresponding fitting windows
	## and fitting windows without corresponding time series
	lts <- intersect(levels(frame_ts$ts), levels(frame_windows$ts))
	if (length(lts) == 0L)
		stop(gettextf("found %d fitting windows with corresponding time series data",
		              0L),
		     domain = NA)
	frame_ts$ts <- factor(frame_ts$ts, levels = lts, exclude = NULL)
	i1 <- !is.na(frame_ts$ts)
	frame_ts <- frame_ts[i1, , drop = FALSE]
	frame_windows$ts <- factor(frame_windows$ts, levels = lts, exclude = NULL)
	i2 <- !is.na(frame_windows$ts)
	frame_windows <- frame_windows[i2, , drop = FALSE]
	frame_parameters <- lapply(frame_parameters, `[`, i2, , drop = FALSE)
	frame_extra <- frame_extra[i2, , drop = FALSE]


	### Another validation step ################################################

	### Time series stuff

	if (!all(is.finite(frame_ts$time)))
		stop(gettextf("'%s' is not finite", nf1$time),
		     domain = NA)
	check_ok_diff_time <-
	if (do_day_of_week <- (model$day_of_week > 0L)) {
		if (is.double(frame_ts$time)) {
			if (!isTRUE(all.equal(frame_ts$time, z <- round(frame_ts$time))))
				stop(gettextf("'%s' is not integer-valued", nf1$time),
				     domain = NA)
			frame_ts$time <- z
		}
		function(x) all(diff(x) == 1)
	}
	else function(x) all(diff(x) > 0)
	if (!all(tapply(frame_ts$time, frame_ts$ts, check_ok_diff_time)))
		stop(switch(1L + 2L * do_day_of_week + nzchar(nf1$ts),
		            gettextf("'%s' must be increasing", nf1$time),
		            gettextf("'%s' must be increasing in each level of %s", nf1$time, nf1$ts),
		            gettextf("'%s' must be increasing with unit spacing", nf1$time, nf1$time),
		            gettextf("'%s' must be increasing with unit spacing in each level of %s", nf1$time, nf1$ts)),
		     domain = NA)
	if (!all(frame_ts$x[!is.na(frame_ts$x)] >= 0))
		stop(gettextf("'%s' is not non-negative", nf1$x),
		     domain = NA)
	if (is.double(frame_ts$x)) {
		is_NaN_or_Inf <- is.nan(frame_ts$x) | is.infinite(frame_ts$x)
		if (any(is_NaN_or_Inf)) {
			warning(gettextf("replacing %f, %f, and %f in '%s' with %d",
			                 NaN, Inf, -Inf, nf1$x, NA),
			        domain = NA)
			frame_ts$x[is_NaN_or_Inf] <- NA
		}
		if (!isTRUE(all.equal(frame_ts$x, z <- round(frame_ts$x))))
			warning(gettextf("replacing '%s' with %s(%s)",
			                 nf1$x, "round", nf1$x),
			        domain = NA)
		frame_ts$x <- z
	}

	### Mixed effects stuff

	check_ok_bar_lhs_double <-
	function(data, names) {
		data <- data[names]
		i <- vapply(data, is.double, FALSE)
		f <- function(x) any(is.infinite(x))
		!any(vapply(data[i], f, FALSE))
	}
	if (!all(mapply(check_ok_bar_lhs_double, data = frame_parameters, names = names_bar_lhs)))
		stop(gettextf("numeric variables on left hand side of '%f' in '%f' contain %f or %f",
		              "|", "formula_parameters", Inf, -Inf),
             domain = NA)


	### Labeling step ##########################################################

	## Order everything by time series and chronologically within time series
	o1 <- order(frame_ts$ts)
	frame_ts <- frame_ts[o1, , drop = FALSE]
	o2 <- do.call(order, unname(frame_windows))
	frame_windows <- frame_windows[o2, , drop = FALSE]
	frame_parameters <- lapply(frame_parameters, `[`, o2, , drop = FALSE)
	frame_extra <- frame_extra[o2, , drop = FALSE]

	## Create enumerated labels for fitting windows as ordered
	N <- nrow(frame_windows)
	lw <- sprintf("window_%0*d", nchar(sprintf("%d", N)), seq_len(N))
	frame_windows$window <- gl(N, 1L, labels = lw)

	## Create a factor grouping observations by fitting window
	make_window_segment <-
	function(d1, d2) {

		if (!all(d2$start < d2$end))
			stop(switch(1L + nzchar(nf1$ts),
			            gettextf("fitting windows (%s, %s] in time series \"%s\" do not satisfy %s < %s",
			                     nf2$start, nf2$end, d2$ts[1L], nf2$start, nf2$end),
			            gettextf("fitting windows (%s, %s] do not satisfy %s < %s",
			                     nf2$start, nf2$end,            nf2$start, nf2$end)),
			     domain = NA)
		if (!all(d2$start[-1L] >= d2$end[-nrow(d2)]))
			stop(switch(1L + nzchar(nf1$ts),
			            gettextf("fitting windows (%s, %s] in time series \"%s\" are not disjoint",
			                     nf2$start, nf2$end, d2$ts[1L]),
			            gettextf("fitting windows (%s, %s] are not disjoint",
			                     nf2$start, nf2$end           )),
			     domain = NA)
		f <- function(a, b) which(d1$time >= a & d1$time <= b)[-1L]
		index <- Map(f, a = d2$start, b = d2$end)
		ulindex <- unlist1(index)
		if (na_action_ts == "fail" && anyNA(d1$x[ulindex]))
			stop(switch(1L + nzchar(nf1$ts),
			            gettextf("fitting windows (%s, %s] in time series \"%s\" contain %d in '%s'",
			                     nf2$start, nf2$end, d2$ts[1L], NA, nf1$x),
			            gettextf("fitting windows (%s, %s] contain %d in '%s'",
			                     nf2$start, nf2$end,            NA, nf1$x)),
			     domain = NA)
		if (any(vapply(index, function(i) sum(!is.na(d1$x[i])) == 0L, FALSE)))
			stop(switch(1L + nzchar(nf1$ts),
			            gettextf("fitting windows (%s, %s] in time series \"%s\" contain %d observations of '%s'",
			                     nf2$start, nf2$end, d2$ts[1L], 0L, nf1$x),
			            gettextf("fitting windows (%s, %s] contain %d observations of '%s'",
			                     nf2$start, nf2$end,            0L, nf1$x)),
			     domain = NA)
		window <- rep.int(factor(NA, levels = levels(d2$window)), nrow(d1))
		window[ulindex] <- rep.int(d2$window, lengths(index))
		window
	}
	window_split <- Map(make_window_segment,
	                    d1 = split(frame_ts, frame_ts$ts),
	                    d2 = split(frame_windows, frame_windows$ts))
	frame_ts$window <- unsplit(window_split, frame_ts$ts)


	### Cleaning up ############################################################

	## Ignore edge cases
	frame_ts$x[!duplicated(frame_ts$ts)] <- NA

	## Contract fitting windows to narrowest intervals
	## containing the same set of observations
	first <- which(!(is.na(frame_ts$window) | duplicated(frame_ts$window))) - 1L
	last  <- which(!(is.na(frame_ts$window) | duplicated(frame_ts$window, fromLast = TRUE)))
	frame_windows$start <- frame_ts$time[first]
	frame_windows$end   <- frame_ts$time[last]

	frame_ts <- frame_ts[c("ts", "window", "time", "x")]
	attr(frame_ts, "first") <- first
	attr(frame_ts, "last") <- last
	row.names(frame_ts) <- NULL

	frame_windows <- frame_windows[c("ts", "window", "start", "end")]
	row.names(frame_windows) <- NULL

	for (i in seq_along(frame_parameters)) {
		frame_parameters[[i]] <- droplevels(frame_parameters[[i]])
		attr(frame_parameters[[i]], "terms") <- terms(formula_parameters[[i]])
		row.names(frame_parameters[[i]]) <- NULL
	}

	frame_extra <- droplevels(frame_extra)
	row.names(frame_extra) <- NULL

	list(ts = frame_ts,
	     windows = frame_windows,
	     parameters = frame_parameters,
	     extra = frame_extra)
}

egf_make_priors <-
function(formula_priors, top, beta, theta, Sigma) {
	l <- list(beta = beta, theta = theta, Sigma = Sigma)
	len <- vapply(l, `[[`, 0L, "length")
	index_full <- lapply(len, seq_len)

	n <- length(formula_priors)
	p <- length(top$names)
	nms <- names(l)[len > 0L]

	priors <- list(top = vector("list", p),
	               bottom = lapply(len, vector, mode = "list"))
	names(priors$top) <- top$names

	for (i in seq_len(n)) {
		x <- formula_priors[[i]]

		lhs <- x[[2L]]
		if (is.name(lhs) && any(m <- (lhs == nms))) {
			kind <- nms[m]
			index <- index_full[[nms[m]]]
		}
		else if (is.call(lhs) &&
		         (lhs[[1L]] == "[" || lhs[[1L]] == "[[") &&
		         any(m <- (lhs[[2L]] == nms))) {
			kind <- nms[m]
			index <- eval(lhs, index_full, environment(x))
			if (anyNA(index))
				stop(gettextf("argument of '%s' on left hand side of %s[[%d]] not a valid index vector for '%s' of length %d",
				              deparse(lhs[[1L]]), "formula_priors", i, kind, len[[kind]]),
				     domain = NA)
		}
		else if (any(m <- (lhs == names(priors$top))))
			kind <- names(priors$top)[m]
		else
			stop(gettextf("left hand side of %s[[%d]] does not match any accepted format",
			              "formula_priors", i),
			     domain = NA)

		rhs <- x[[3L]]
		prior <- eval(rhs, environment(x))
		if (!inherits(prior, "egf_prior"))
			stop(gettextf("right hand side of %s[[%d]] does not evaluate to an \"%s\" object",
			              "formula_priors", i, "egf_prior"),
			     domain = NA)
		if (kind %in% nms) {
			allowed_prior_families <- l[[kind]]$family
			prior$parameters <- lapply(prior$parameters, rep_len,
			                           length.out = length(index))
			f <-
			function(k) {
				prior$parameters <- lapply(prior$parameters, `[[`, k)
				prior
			}
			priors$bottom[[kind]][index] <- lapply(seq_along(index), f)
		}
		else {
			allowed_prior_families <- top$family
			prior$parameters <- lapply(prior$parameters, `[[`, 1L)
			priors$top[[kind]] <- prior
		}
		if (!(prior$family %in% allowed_prior_families))
			stop(gettextf("priors on '%d' are restricted to families %s",
			              kind, deparse(allowed_prior_families)),
			     domain = NA)
	}

	x <- priors$bottom$Sigma
	for (i in seq_along(x)) {
		if (is.null(x[[i]]) || !grepl("^(inv)?wishart$", x[[i]]$family))
			next
		n1 <- as.integer(choose(length(x[[i]]$parameters$scale) + 1L, 2L))
		n2 <- Sigma$rows[i]
		if (n1 != n2)
			stop(gettextf("prior on %s[[%d]] specifies '%s' having %d rows, but the expected number is %d",
			              "Sigma", i, "scale", n1, n2),
			     domain = NA)
	}

	priors
}

##' Construct design matrices
##'
##' Utilities for constructing the design matrix \code{X} and \code{Z}
##' associated with a fixed effects model formula \code{~tt} or random
##' effects term \code{(tt | g)}, respectively.
##'
##' @param fixed
##'   A fixed effects formula \code{~tt}.
##' @param random
##'   A random effects term \code{(tt | g)}.
##'   More precisely, a \link{call} to binary operator \code{`|`}.
##' @param data
##'   A \link[=model.frame]{model frame} listing the variables
##'   used in \code{fixed} or \code{random}.
##' @param sparse
##'   A logical flag. If \code{TRUE}, then the design matrix
##'   is returned in \link[Matrix:sparseMatrix]{sparse} format.
##'
##' @details
##' \code{egf_make_X(fixed, data, sparse)}
##' constructs an \code{X} matrix by evaluating
##' \code{\link{model.matrix}(fixed, data)} or
##' \code{\link[Matrix]{sparse.model.matrix}(fixed, data)}
##' (depending on \code{sparse})
##' and deleting from the result columns containing only zeros.
##'
##' \code{egf_make_Z(x, data)} constructs a \code{Z} matrix
##' following the steps outlined in \code{\link{vignette}("lmer", "lme4")}.
##' It uses \code{egf_make_X} to construct the so-called raw model matrix
##' from \code{random[[2L]]}.
##' The result is always sparse.
##'
##' @return
##' A \link[=matrix]{dense} or \link[Matrix:sparseMatrix]{sparse} matrix
##' with attributes:
##' \item{assign}{
##'   See \code{\link{model.matrix}}.
##'   Indexes \code{\link{labels}(\link{terms}(~tt))}
##'   for the expression \code{tt} in \code{fixed ~ tt}
##'   or \code{random = (tt | g)}.
##' }
##' \item{contrasts}{
##'   See \code{\link{model.matrix}}.
##'   Absent if the expression \code{tt} in \code{fixed = ~tt}
##'   or \code{random = (tt | g)} does not contain terms
##'   evaluating to character vectors or factors.
##' }
##' \item{level}{
##'   (\code{egf_make_Z} only.)
##'   A factor of length \code{ncol(Z)} with levels \code{levels(g)},
##'   useful for splitting the columns of \code{Z} by group level.
##' }

egf_make_X <-
function(fixed, data, sparse) {
	mm <- if (sparse) sparse.model.matrix else model.matrix
	X <- mm(fixed, data = data)
	j <- colSums(abs(X)) > 0
	structure(X[, j, drop = FALSE],
	          assign = attr(X, "assign")[j],
	          contrasts = attr(X, "contrasts"))
}

egf_make_Z <-
function(random, data) {
	X <- egf_make_X(as.formula(call("~", random[[2L]])), data = data, sparse = FALSE)
	ng <- vapply(split_interaction(random[[3L]]), deparse1, "")
	g <- interaction(data[ng], drop = TRUE, sep = ":", lex.order = FALSE)
	J <- as(g, Class = "sparseMatrix")
	Z <- t(KhatriRao(J, t(X)))
	## For consistency, we desire 'model.matrix'-style names for group levels
	G <- model.matrix(as.formula(call("~", call("+", 0, random[[3L]]))), data = data)
	j <- colSums(abs(G)) > 0
	G <- G[, j, drop = FALSE]
	rownames(Z) <- as.character(seq_len(nrow(Z)))
	colnames(Z) <- sprintf("(%s | %s)",
	                       rep(colnames(X), times = ncol(G)),
	                       rep(colnames(G), each = ncol(X)))
	structure(Z,
	          assign = rep(attr(X, "assign"), times = ncol(G)),
	          contrasts = attr(X, "contrasts"),
	          level = gl(ncol(G), ncol(X), labels = levels(g)))
}

##' Combine design matrices
##'
##' Utilities for combining parameter-specific fixed effects design matrices
##' \code{X} and term-specific random effects design matrices \code{Z}.
##'
##' @param fixed
##'   A named list of fixed effects formulae \code{~tt}.
##'   \code{names(fixed)} must indicate corresponding top level
##'   nonlinear model parameters.
##' @param X
##'   A list of \code{X} matrices obtained by applying
##'   \code{egf_make_X} to the elements of \code{fixed}.
##' @param random
##'   A named list of random effects terms \code{(tt | g)}.
##'   More precisely, a named list of \link{call}s to binary operator
##'   \code{`|`}.
##'   \code{names(random)} must indicate corresponding top level
##'   nonlinear model parameters.
##' @param Z
##'   A list of \code{Z} matrices obtained by applying \code{egf_make_Z}
##'   to the elements of \code{random}.
##'
##' @return
##' \code{egf_combine_X} returns a list with elements:
##' \item{X}{
##'   The result of combining (in the sense of \code{\link{cbind}})
##'   the supplied design matrices.
##' }
##' \item{effects}{
##'   A data frame with one row per column of the combined design matrix,
##'   storing details about the corresponding linear coefficients.
##' }
##' \item{contrasts}{
##'   The result of concatenating the so-named attributes of the supplied
##'   design matrices and deleting duplicates.
##' }
##' \code{egf_combine_Z} returns a similar list, with \code{Z} replacing
##' \code{X}. However, the columns of \code{Z} and rows of \code{effects}
##' are permuted so that they are ordered
##' by relation to a common random effects term \code{(tt | g)}
##' (order of terms is by appearance in \code{random}), then
##' by relation to a common level of grouping variable \code{g}
##' (order of levels of interactions is reverse lexicographic), then
##' by top level nonlinear model parameter
##' (order of parameters is by appearance in \code{names(random)}).
##'
##' Element \code{effects} stores the following details about coefficients:
##' \item{cov}{
##'   (\code{egf_combine_Z} only.)
##'   Name of a covariance matrix, in the format \code{"Sigma[\%d]"}.
##'   This is the interaction of \code{term} and \code{group}.
##' }
##' \item{vec}{
##'   (\code{egf_combine_Z} only.)
##'   Name of a random vector, in the format \code{"u[\%d]"}
##'   This is the interaction of \code{term}, \code{group}, and \code{level}.
##' }
##' \item{bottom}{
##'   Name of a bottom level mixed effects model parameter,
##'   in the format \code{"beta[\%d]"} or \code{"b[\%d]"}.
##' }
##' \item{top}{
##'   Name of the top level nonlinear model parameter whose fitted value
##'   is a function of \code{bottom}.
##' }
##' \item{term}{
##'   Deparsed term from \code{tt}, or \code{"(Intercept)"}.
##' }
##' \item{group}{
##'   (\code{egf_combine_Z} only.)
##'   Deparsed expression \code{g}.
##' }
##' \item{level}{
##'   (\code{egf_combine_Z} only.)
##'   Level of factor or interaction indicated by \code{group}.
##' }
##' \item{colname}{
##'   Column name in design matrix.
##' }

egf_combine_X <-
function(fixed, X) {
	stopifnot(length(fixed) == length(X))

	nc <- vapply(X, ncol, 0L)
	assign <- lapply(X, attr, "assign")
	contrasts <- lapply(X, attr, "contrasts")

	get_term_labels <-
	function(formula, assign)
		c("(Intercept)", labels(terms(formula)))[1L + assign]

	bottom <- disambiguate(rep.int("beta", sum(nc)))
	top <- rep.int(gl(length(fixed), 1L, labels = names(fixed)), nc)
	term <- Map(get_term_labels, formula = fixed, assign = assign)
	term <- factor(unlist1(term))

	X <- do.call(cbind, X)
	effects <- data.frame(bottom, top, term, colname = colnames(X),
	                      row.names = NULL, stringsAsFactors = FALSE)
	contrasts <- unlist(unname(contrasts), recursive = FALSE, use.names = TRUE)
	contrasts[duplicated(names(contrasts))] <- NULL
	list(X = X, effects = effects, contrasts = contrasts)
}

egf_combine_Z <-
function(random, Z) {
	stopifnot(length(random) == length(Z))
	if (length(random) == 0L)
		return(NULL)

	random_formula <- lapply(random, function(x) as.formula(call("~", x[[2L]])))
	random_group <- lapply(random, `[[`, 3L)

	nc <- vapply(Z, ncol, 0L)
	assign <- lapply(Z, attr, "assign")
	contrasts <- lapply(Z, attr, "contrasts")
	level <- lapply(Z, attr, "level")

	get_term_labels <-
	function(formula, assign)
		c("(Intercept)", labels(terms(formula)))[1L + assign]

	bottom <- disambiguate(rep.int("b", sum(nc)))
	top <- factor(rep.int(names(random), nc), levels = unique(names(random)))
	term <- Map(get_term_labels, formula = random_formula, assign = assign)
	term <- factor(unlist1(term))
	group <- factor(rep.int(vapply(random_group, deparse1, ""), nc))
	level <- unlist1(level)
	cov <- interaction(term, group, drop = TRUE, lex.order = TRUE)
	levels(cov) <- disambiguate(rep.int("Sigma", nlevels(cov)))
	vec <- interaction(term, group, level, drop = TRUE, lex.order = TRUE)
	levels(vec) <- disambiguate(rep.int("u", nlevels(vec)))

	Z <- do.call(cbind, Z)
	effects <- data.frame(cov, vec, bottom, top, term, group, level,
	                      colname = colnames(Z),
	                      row.names = NULL, stringsAsFactors = FALSE)
	contrasts <- unlist(unname(contrasts), recursive = FALSE, use.names = TRUE)
	contrasts[duplicated(names(contrasts))] <- NULL

	o <- do.call(order, unname(effects[c("cov", "vec", "top")]))
	Z <- Z[, o, drop = FALSE]
	effects[-3L] <- effects[o, -3L, drop = FALSE]
	list(Z = Z, effects = effects, contrasts = contrasts)
}
