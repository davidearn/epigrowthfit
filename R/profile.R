profile.egf <-
function(fitted,
         A,
         which,
         level = 0.95,
         grid = 12L,
         parallel = egf_parallel(),
         trace = FALSE,
         top = egf_top(fitted),
         subset = NULL,
         select = NULL,
         ...) {
	stopifnot(sum(m. <- c(missing(A), missing(which), missing(top))) >= 2L,
	          isNumber(level), level > 0, level < 1,
	          isNumber(grid), grid >= 1,
	          inherits(parallel, "egf_parallel"),
	          isTrueFalse(trace))

	par <- coef(fitted)

	## If profiling user-specified linear combinations
	## of elements of c(beta, theta)
	if (!is.null(A)) {
		if (is.double(A) && !is.object(A) && !is.array(A))
			dim(A) <- c(1L, length(A))
		stopifnot((is.double(A) && !is.object(A) && is.matrix(A)) ||
		          is(A, "dMatrix"),
		          nrow(A) > 0L,
		          ncol(A) == length(par),
		          is.finite(range(A)),
		          min(rowSums(A != 0)) > 0)
	}
	## If profiling user-specified elements of c(beta, theta)
	else if (!is.null(which)) {
		which <- `names<-`(seq_along(par), labels(par))[which]
		A <- new("dgRMatrix",
	             Dim = c(length(which), length(par)),
	             Dimnames = list(names(which), NULL),
	             p = 0:length(which),
	             j = which - 1L,
	             x = rep.int(1, length(which)))
	}
	## If profiling population fitted values
	## of top level nonlinear model parameters
	else if (!is.null(top)) {
		top. <- egf_top(fitted)
		top <- unique(match.arg(top, top., several.ok = TRUE))

		frame_windows <- model.frame(fitted, "windows")
		frame_combined <- model.frame(fitted, "combined")
		subset <- egf_eval_subset(subset, frame_combined, parent.frame())
		select <- egf_eval_select(select, frame_combined, baseenv())

		l <- egf_preprofile(fitted, subset = subset, top = top)
		Y <- l[["Y"]]
		A <- l[["A"]]
	}
	else stop(gettextf("one of '%s', '%s', and '%s' must be non-NULL",
	                   "A", "which", "top"),
	          domain = NA)

	## Covariance matrix of c(beta, theta)
	V <- unclass(vcov(fitted))
	m <- nrow(A)
	if (!m.[2L]) {
		a <- which
		h <- 0.25 * sqrt(diag(V, names = FALSE)[which])
		s <- "name"
	}
	else {
		## Covariance matrix of A %*% c(beta, theta)
		V <- A %*% tcrossprod(V, A)
		a <- lapply(seq_len(m), function(i) A[i, ])
		h <- 0.25 * sqrt(diag(V, names = FALSE))
		s <- "lincomb"
	}

	## y := nll_restricted - nll_minimum = 0.5 * (change in deviance)
	ytol <- 0.5 * qchisq(level, df = 1)
	ystep <- ytol / grid
	nomp <- fitted[["control"]][["omp_num_threads"]]

	doProfile <-
	function(i, a, h) {
		## FIXME: write own version of 'tmbprofile' giving proper 'par.vals'!
		if (trace)
			cat(gettextf("computing likelihood profile %d of %d ...",
			             i, m),
			    "\n", sep = "")
		args <- list(obj = obj, h = h, ytol = ytol, ystep = ystep)
		args[[s]] <- a
		ans <- do.call(TMB::tmbprofile, args)
		dot <- ans[[1L]] # it's a "dot" product
		nll <- ans[[2L]]
		i.min <- which.min(nll) # this row is duplicated => [-i.min] below
		z <- sqrt(2 * (nll - nll[i.min]))
		z[seq_len(i.min)] <- -z[seq_len(i.min)]
		ans <- list(z[-i.min], as.matrix(dot[-i.min]))
		attributes(ans) <- list(names = c("z", "par.vals"),
		                        row.names = .set_row_names(length(z) - 1L),
		                        class = "data.frame")
		ans
	}

	if (parallel[["method"]] == "snow") {
		environment(doProfile) <- globalenv()

		## Reconstruct list of arguments to 'MakeADFun' from object internals
		## for retaping
		args <- egf_tmb_remake_args(fitted[["tmb_out"]], par = fitted[["best"]])

		## Retrieve path to shared object for loading
		dll <- .dll

		cl <- parallel[["cl"]]
		if (is.null(cl)) {
			cl <- do.call(makePSOCKcluster, parallel[["args"]])
			on.exit(stopCluster(cl), add = TRUE)
		}
		vars <- c("dll", "nomp", "args", "trace", "m", "s", "ytol", "ystep")
		clusterExport(cl, varlist = vars, envir = environment())
		clusterEvalQ(cl, {
			dyn.load(dll)
			if (TMB::openmp(n = NULL) > 0L)
				TMB::openmp(n = nomp)
			obj <- do.call(TMB::MakeADFun, args)
		})
		ans <- clusterMap(cl, doProfile, i = seq_len(m), a = a, h = h)
	}
	else {
		if (nzchar(parallel[["outfile"]])) {
			outfile <- file(parallel[["outfile"]], open = "wt")
			sink(outfile, type = "output")
			sink(outfile, type = "message")
			on.exit(add = TRUE, {
				sink(type = "message")
				sink(type = "output")
			})
		}
		onomp <- TMB::openmp(n = NULL)
		if (onomp > 0L) {
			TMB::openmp(n = nomp)
			on.exit(TMB::openmp(n = onomp), add = TRUE)
		}
		obj <- fitted[["tmb_out"]]
		ans <- switch(parallel[["method"]],
		              multicore = do.call(mcMap, c(list(f = doProfile, i = seq_len(m), a = a, h = h), parallel[["args"]])),
		              serial = Map(doProfile, i = seq_len(m), a = a, h = h))
	}

	if (m.[1L] && m.[2L]) {
		attr(ans, "frame") <-
			data.frame(frame_windows[subset, c("ts", "window"), drop = FALSE],
			           frame_combined[subset, select, drop = FALSE],
			           row.names = NULL,
			           check.names = FALSE,
			           stringsAsFactors = FALSE)
		attr(ans, "top") <- factor(top, levels = top.)
		attr(ans, "Y") <- Y
	}
	names(ans) <-
		if (is.null(rn <- dimnames(A)[[1L]]))
			paste0("(A %%*%% coef(.))", "[", format(seq_along(ans)), "]")
		else rn
	attr(ans, "A") <- A
	attr(ans, "level") <- level
	class(ans) <- c("profile.egf", "profile")
	ans
}

confint.profile.egf <-
function(object, parm = seq_along(object), level = attr(object, "level"),
         maybeClass = FALSE, ...) {
	stopifnot(isNumber(level), level > 0, level < 1,
	          level <= attr(object, "level"))
	which <- `names<-`(seq_along(object), names(object))[parm]
	if (anyNA(which))
		stop(gettextf("invalid '%s'", "parm"), domain = NA)
	h <- 0.5 * (1 - level)
	p <- c(h, 1 - h)
	q <- qnorm(p)
	percent <- formatp(p)
	if (doClass <- maybeClass) {
		l <- attributes(object)[c("frame", "top", "Y")]
		doClass <- !any(vapply(l, is.null, FALSE))
	}
	ans <- array(NA_real_,
	             dim = c(length(which), 2L),
	             dimnames = list(if (!doClass) names(which), percent))
	for (i in which) {
		pr <- object[[i]]
		sp <- spline(x = pr[["par.vals"]][, 1L], y = pr[["z"]])
		ans[i, ] <- approx(x = sp[["y"]], y = sp[["x"]], xout = q)[["y"]]
	}
	if (doClass) {
		dY <- dim(l[["Y"]])
		ans <- data.frame(
			top = rep(l[["top"]], each = dY[1L])[which],
			l[["frame"]][rep.int(seq_len(dY[1L]), dY[2L])[which], ],
			ci = I(ans + l[["Y"]][which]),
			row.names = NULL,
			check.names = FALSE,
			stringsAsFactors = FALSE)
		class(ans) <- c("confint.egf", oldClass(ans))
	}
	ans
}

plot.profile.egf <-
function(x, parm = seq_along(x), level = attr(x, "level"),
         type = c("z", "abs(z)", "z^2"), ...) {
	type <- match.arg(type)
	if (!is.null(level))
	stopifnot(is.double(level), !anyNA(level),
	          min(1, level) > 0, max(0, level) < 1,
	          max(0, level) <= attr(x, "level"))
	which <- `names<-`(seq_along(x), names(x) -> nms)[parm]
	if (anyNA(which))
		stop(gettextf("invalid '%s'", "parm"), domain = NA)
	h <- 0.5 * (1 - attr(x, "level"))
	p <- c(h, 1 - h)
	q <- qnorm(p)
	ylim <- range(vapply(x, function(pr) range(pr[["z"]]), double(2L)))
	ylim <- switch(type,
	               "z"      = ylim,
	               "abs(z)" = c(0, ylim[2L]),
	               "z^2"    = c(0, ylim[2L]^2))
	ylab <- switch(type,
	               "z"      = quote(z),
	               "abs(z)" = quote(group("|", z, "|")),
	               "z^2"    = quote(z^2))
	h <- 0.5 * (1 - level)
	p <- rbind(h, 1 - h, deparse.level = 0L)
	q <- qnorm(p)
	for (i in which) {
		dev.hold()
		pr <- x[[i]]
		pr.x <- pr[["par.vals"]][, 1L]
		pr.y <- pr[["z"]]
		pr.y. <- switch(type,
		                "z"      = ,
		                "abs(z)" = pr.y,
		                "z^2"    = sign(pr.y) * pr.y^2)
		sp <- spline(x = pr.x, y = pr.y.)
		sp.x <- sp[["x"]]
		sp.y <- sp[["y"]]
		sp.y. <- switch(type,
		                "z"      = sp.y,
		                "abs(z)" = ,
		                "z^2"    = abs(sp.y))
		plot(x = sp.x, y = sp.y., type = "n",
		     xlim = NULL, ylim = ylim, xlab = nms[i], ylab = ylab, ...)
		lines(x = sp.x, y = sp.y., lty = 1L)
		usr <- par("usr")
		for (j in seq_along(level)) {
			pj <- p[, j]
			qj <- q[, j]
			qj <- switch(type,
			             "z"      = ,
			             "abs(z)" = qj,
			             "z^2"    = sign(qj) * qj^2)
			hj <- switch(type,
			             "z"      = qj,
			             "abs(z)" = ,
			             "z^2"    = c(usr[3L], qj[2L]))
			vj <- approx(x = sp.y., y = sp.x, xout = qj)[["y"]]
			segments(x0 = vj[c(1L, 1L, 2L, if (type == "z") 2L)],
			         x1 = vj[c(1L, 2L, 2L, if (type == "z") 1L)],
			         y0 = hj[c(1L, 2L, 2L, if (type == "z") 1L)],
			         y1 = hj[c(2L, 2L, 1L, if (type == "z") 1L)],
			         lty = 3L)
			if (type == "z") {
			abline(h = 0, lty = 3L)
			text(x = vj[1L],
			     y = hj[2L],
			     labels = formatp(pj[2L]),
			     adj = c(0, 0))
			text(x = vj[2L],
			     y = hj[1L],
			     labels = formatp(pj[1L]),
			     adj = c(1, 0))
			}
			else
			text(x = mean(usr[1:2]),
			     y = hj[2L],
			     labels = formatp(pj[2L] - pj[1L]),
			     adj = c(0.5, 0))
		}
		dev.flush()
	}
	invisible(NULL)
}
