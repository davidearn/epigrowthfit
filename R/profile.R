profile.egf <-
function(fitted,
         level = 0.95,
         A = seq_along(par),
         grid = 12L,
         parallel = egf_parallel(),
         trace = FALSE,
         top = egf_top(fitted),
         subset = NULL,
         select = NULL,
         ...) {
	stopifnot(isNumber(level), level > 0, level < 1,
	          isNumber(grid), grid >= 1,
	          inherits(parallel, "egf_parallel"),
	          isTrueFalse(trace))

	par <- coef(fitted)
	which <- NULL

	if (m.A <- is.null(A)) {
		top. <- egf_top(fitted)
		top <- unique(match.arg(top, top., several.ok = TRUE))

		frame <- model.frame(fitted, "combined")
		subset <- egf_eval_subset(subset, frame, parent.frame())
		select <- egf_eval_select(select, frame, baseenv())
		frame <- frame[subset, select]

		info <- model.frame(fitted, "windows")[subset, c("ts", "window")]
		ts     <- info[["ts"    ]]
		window <- info[["window"]]

		A <- egf_make_A(fitted, top = top, subset = subset)
	}
	else if ((is.double(A) && !is.object(A) && is.matrix(A)) ||
	         is(A, "dMatrix"))
		stopifnot(nrow(A) > 0L,
		          ncol(A) == 1L + length(par),
		          is.finite(range(A)),
		          min(rowSums(A[, -1L, drop = FALSE] != 0)) > 0)
	else {
		which <- `names<-`(seq_along(par), labels(par))[A]
		A <- new("dgRMatrix",
		         Dim = c(length(which), 1L + length(par)),
		         Dimnames = list(names(which), NULL),
		         p = c(0L, seq_along(which)),
		         j = which,
		         x = rep.int(1, length(which)))
	}
	a <- lapply(seq_len(nrow(A)), function(i) A[i, ])

	V. <- vcov(fitted)
	if (!is.null(which))
		h <- 0.25 * sqrt(diag(V., names = FALSE)[which])
	else {
		A. <- A[, -1L, drop = FALSE]
		V. <- A. %*% tcrossprod(V., A.)
		h <- 0.25 * sqrt(diag(V., names = FALSE))
	}

	## y := nll_restricted - nll_minimum = 0.5 * (change in deviance)
	ytol <- 0.5 * qchisq(level, df = 1)
	ystep <- ytol / grid

	nomp <- fitted[["control"]][["omp_num_threads"]]
	nprf <- length(a)

	doProfile <-
	function(i, a, h) {
		## FIXME: write own version of 'tmbprofile' giving proper 'par.vals'!
		if (trace)
			cat(gettextf("computing profile %d of %d ...", i, nprf),
			    "\n", sep = "")
		ans <- TMB::tmbprofile(obj = obj, lincomb = a[-1L], h = h,
		                       ytol = ytol, ystep = ystep)
		dot <- ans[[1L]] # it's a "dot" product
		nll <- ans[[2L]]
		i.min <- which.min(nll) # this row is duplicated => [-i.min] below
		z <- sqrt(2 * (nll - nll[i.min]))
		z[seq_len(i.min)] <- -z[seq_len(i.min)]
		ans <- list(z[-i.min],
		            as.matrix((if (is.integer(a)) 0 else a[1L]) + dot[-i.min]))
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
		vars <- c("dll", "nomp", "nprf", "args", "trace",
		          "ytol", "ystep")
		clusterExport(cl, varlist = vars, envir = environment())
		clusterEvalQ(cl, {
			dyn.load(dll)
			if (TMB::openmp(n = NULL) > 0L)
				TMB::openmp(n = nomp)
			obj <- do.call(TMB::MakeADFun, args)
		})
		ans <- clusterMap(cl, doProfile, i = seq_len(nprf), a = a, h = h)
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
		              multicore = do.call(mcMap, c(list(f = doProfile, i = seq_len(nprf), a = a, h = h), parallel[["args"]])),
		              serial = Map(doProfile, i = seq_len(nprf), a = a, h = h))
	}

	if (m.A) {
		dim(ans) <- c(length(subset), length(top))
		dimnames(ans) <- list(rn <- paste(ts, window, sep = ", "), top)
		names(ans) <- paste(rep(top, each = length(subset)), rn, sep = ", ")
		attr(ans, "top") <- factor(top, levels = top.)
		attr(ans, "ts") <- ts
		attr(ans, "window") <- window
		attr(ans, "frame") <- frame
	}
	else if (is.null(rn <- dimnames(A)[[1L]]))
		names(ans) <- paste0("(A %%*%% c(1, coef(.)))",
		                     "[", format(seq_along(ans)), "]")
	else
		names(ans) <- rn
	attr(ans, "A") <- A
	attr(ans, "par") <- par
	attr(ans, "level") <- level
	class(ans) <- c("profile.egf", "profile")
	ans
}

confint.profile.egf <-
function(object, parm = seq_along(object), level = attr(object, "level"),
         class = FALSE, ...) {
	stopifnot(isNumber(level), level > 0, level < 1,
	          level <= attr(object, "level"),
	          isTrueFalse(class))
	m.A <- !is.null(attr(object, "frame"))
	which <- `names<-`(seq_along(object), names(object))[parm]
	if (anyNA(which))
		stop(gettextf("invalid '%s'", "parm"), domain = NA)

	h <- 0.5 * (1 - level)
	p <- c(h, 1 - h)
	q <- qnorm(p)
	percent <- formatp(p, 3L)
	if (class)
		class <- m.A

	ans <- array(NA_real_, dim = c(length(which), 2L))
	for (i in which) {
		pr <- object[[i]]
		sp <- spline(x = pr[["par.vals"]][, 1L], y = pr[["z"]])
		ans[i, ] <- approx(x = sp[["y"]], y = sp[["x"]], xout = q)[["y"]]
	}

	if (class || !m.A || !missing(parm))
		dimnames(ans) <- list(if (!class) names(which), percent)
	else {
		dim(ans) <- c(dim(object), 2L)
		dimnames(ans) <- c(dimnames(object), list(percent))
	}

	if (class) {
		a. <- attributes(object)
		ns <- nrow(a.[["frame"]])
		nt <- length(a.[["top"]])
		i. <- rep.int(seq_len(ns), nt)[which]
		j. <- rep(seq_len(nt), each = ns)[which]
		tmp <- data.frame(top = a.[["top"]][j.],
		                  ts = a.[["ts"]][i.],
		                  window = a.[["window"]][i.],
		                  value = as.vector(a.[["A"]] %*% c(1, a.[["par"]])),
		                  ci = NA,
		                  a.[["frame"]][i., ],
		                  row.names = NULL,
		                  check.names = FALSE,
		                  stringsAsFactors = FALSE)
		tmp[["ci"]] <- ans
		ans <- tmp
		attr(ans, "level") <- level
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
			vj <- approx(x = sp.y, y = sp.x, xout = qj)[["y"]]
			segments(x0 = vj[c(1L, 1L, 2L, if (type == "z") 2L)],
			         x1 = vj[c(1L, 2L, 2L, if (type == "z") 1L)],
			         y0 = hj[c(1L, 2L, 2L, if (type == "z") 1L)],
			         y1 = hj[c(2L, 2L, 1L, if (type == "z") 1L)],
			         lty = 3L)
			if (type == "z") {
			abline(h = 0, lty = 3L)
			text(x = vj[1L],
			     y = hj[2L],
			     labels = formatp(pj[2L], 3L),
			     adj = c(0, 0))
			text(x = vj[2L],
			     y = hj[1L],
			     labels = formatp(pj[1L], 3L),
			     adj = c(1, 0))
			}
			else
			text(x = mean(usr[1:2]),
			     y = hj[2L],
			     labels = formatp(pj[2L] - pj[1L], 3L),
			     adj = c(0.5, 0))
		}
		dev.flush()
	}
	invisible(NULL)
}
