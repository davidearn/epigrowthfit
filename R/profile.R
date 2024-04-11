profile.egf <-
function(fitted,
         level = 0.95,
         A = NULL,
         which = NULL,
         top = egf_top(fitted),
         grid = 12L,
         parallel = egf_parallel(),
         trace = FALSE,
         subset = NULL,
         select = NULL,
         ...) {
	stopifnot(sum(m. <- c(missing(A), missing(which), missing(top))) >= 2L,
	          is_number_in_interval(level, 0, 1, "()"),
	          is_number_in_interval(grid, 1, Inf, "[)"),
	          inherits(parallel, "egf_parallel"),
	          is_true_or_false(trace))
	n <- sum(!fitted[["random"]])


	## If profiling user-specified linear combinations
	## of elements of c(beta, theta)
	if (!is.null(A)) {
		if (is.double(A) && !is.object(A) && !is.array(A))
			dim(A) <- c(1L, length(A))
		stopifnot((is.double(A) && !is.object(A) && is.matrix(A)) ||
		          is(A, "dMatrix"),
		          nrow(A) > 0L,
		          ncol(A) == n,
		          is.finite(range(A)),
		          min(rowSums(A != 0)) > 0)
	}
	## If profiling user-specified elements of c(beta, theta)
	else if (!is.null(which)) {
		which <- seq_len(n)[which]
		A <- as(new("dgRMatrix",
	                Dim = c(length(which), n),
	                p = 0:(length(which) - 1L),
	                j = which - 1L,
	                x = rep.int(1, length(which))), "CsparseMatrix")
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

	do.profile <-
	function(i, a, h) {
		if (trace)
			cat(gettextf("computing likelihood profile %d of %d ...",
			             i, m),
			    "\n", sep = "")
		args <- list(obj = obj, h = h, ytol = ytol, ystep = ystep)
		args[[s]] <- a
		ans <- do.call(TMB::tmbprofile, args)
		nll <- ans[[2L]]
		i.min <- which.min(nll)
		## (change in deviance) = 2 * (nll_restricted - nll_minimum)
		ans[[2L]] <- 2 * (nll - nll[i.min])
		names(ans) <- c("value", "deviance")
		ans[-i.min, , drop = FALSE] # 'tmbprofile' duplicates this row
	}

	if (parallel[["method"]] == "snow") {
		environment(do.profile) <- globalenv()

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
		ans <- clusterMap(cl, do.profile, i = seq_len(m), a = a, h = h)
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
		              multicore = do.call(mcMap, c(list(f = do.profile, i = seq_len(m), a = a, h = h), parallel[["args"]])),
		              serial = Map(do.profile, i = seq_len(m), a = a, h = h))
	}

	nr <- vapply(ans, nrow, 0L)
	ans <- data.frame(linear_combination = rep.int(gl(m, 1L), nr),
	                  do.call(rbind, ans),
	                  row.names = NULL)
	if (m.[1L] && m.[2L]) {
		i <- rep.int(rep.int(subset, length(top)), nr)
		ans <- data.frame(top = rep.int(rep(factor(top, levels = top.),
		                                    each = length(subset)), nr),
		                  frame_windows[i, c("ts", "window"), drop = FALSE],
		                  ans,
		                  frame_combined[i, select, drop = FALSE],
		                  row.names = NULL,
		                  check.names = FALSE,
		                  stringsAsFactors = FALSE)
		ans[["value"]] <- ans[["value"]] + rep.int(as.double(Y), nr)
	}
	attr(ans, "A") <- A
	attr(ans, "level") <- level
	class(ans) <- c("profile.egf", oldClass(ans))
	ans
}

confint.profile.egf <-
function(object, parm, level = attr(object, "level"), link = TRUE, ...) {
	stopifnot(is_true_or_false(link),
	          is_number_in_interval(level, 0, attr(object, "level"), "(]"))
	q <- qchisq(level, df = 1)

	do.solve <-
	function(d) {
		i.min <- which.min(d[["deviance"]])
		i.left <- seq_len(i.min)
		i.right <- seq.int(i.min, nrow(d))
		estimate <- d[["value"]][i.min]
		lower <- approx(x = d[["deviance"]][i.left],
		                y = d[["value"]][i.left],
		                xout = q)[["y"]]
		upper <- approx(x = d[["deviance"]][i.right],
		                y = d[["value"]][i.right],
		                xout = q)[["y"]]
		d1 <- d[1L, , drop = FALSE]
		m <- match(c("value", "deviance"), names(d1), 0L)
		data.frame(d1[seq_len(m[1L] - 1L)],
		           estimate,
		           lower,
		           upper,
		           d1[-seq_len(m[2L])],
		           row.names = FALSE,
		           check.names = FALSE,
		           stringsAsFactors = FALSE)
	}

	ans <- do.call(rbind,
	               by(object, object[["linear_combination"]], do.solve,
	                  simplify = FALSE))
	if (any(names(object) == "top") && !link) {
		elu <- c("estimate", "lower", "upper")
		f <- lapply(egf_link_extract(levels(ans[["top"]])), egf_link_match,
		            inverse = TRUE)
		ans[elu] <- papply(ans[elu], ans[["top"]], f)
		levels(ans[["top"]]) <- egf_link_remove(levels(ans[["top"]]))
	}
	ans[["linear_combination"]] <- as.integer(as.character(ans[["linear_combination"]]))
	row.names(ans) <- NULL
	attr(ans, "A") <- attr(object, "A")
	attr(ans, "level") <- level
	ans
}

plot.profile.egf <-
function(x, level = attr(x, "level"), sqrt = FALSE, subset = NULL, ...) {
	subset <- egf_eval_subset(subset, x, parent.frame())
	subset <- match(levels(factor(x[["linear_combination"]][subset])),
	                levels(x[["linear_combination"]]))

	stopifnot(is_true_or_false(sqrt))
	f <- if (sqrt) base::sqrt else identity
	do.segments <-
		!sqrt &&
		is.numeric(level) &&
		length(level) > 0L &&
		!all(is.na(level))
	if (do.segments) {
		level <- level[!is.na(level)]
		stopifnot(level > 0, level <= attr(x, "level"))
		## Line segments at heights 'h' in all plots
		h <- qchisq(level, df = 1)
		## Line segment 'j' to start at 'v.lower[[i]][j]'
		## and end at 'v.upper[[i]][j]' in plot 'i'
		ci <- lapply(level, function(p) confint(x, level = p))
		v.lower <- lapply(subset, function(i) vapply(ci, `[`, 0, i, "lower"))
		v.upper <- lapply(subset, function(i) vapply(ci, `[`, 0, i, "upper"))
	}

	dots <- list(...)
	main <-
		if (!is.null(dots[["main"]]))
			rep_len(dots[["main"]], length(subset))
		else if (any(names(x) == "top"))
			c(tapply(as.character(x[["window"]]),
			         x[["linear_combination"]], `[[`, 1L))[subset]
		else character(length(subset))
	xlab <-
		if (!is.null(dots[["xlab"]]))
			rep_len(dots[["xlab"]], length(subset))
		else if (any(names(x) == "top"))
			c(tapply(as.character(x[["top"]]),
			         x[["linear_combination"]], `[[`, 1L))[subset]
		else paste0("linear combination ",
		            levels(x[["linear_combination"]])[subset])
	dots[["main"]] <- dots[["xlab"]] <- NULL
	if (is.null(dots[["ylab"]])) {
		dots[["ylab"]] <- "deviance"
		if (sqrt)
			dots[["ylab"]] <- bquote(sqrt(.(dots[["ylab"]])))
	}
	if (is.null(dots[["ylim"]]))
		dots[["ylim"]] <- c(0, f(max(x[["deviance"]], na.rm = TRUE)))

	for (i in seq_along(subset)) {
		args <- list(formula = f(deviance) ~ value,
		             data = x,
		             subset = (as.integer(x[["linear_combination"]]) == subset[i]),
		             main = main[i],
		             xlab = xlab[i])
		do.call(plot, c(args, dots))
		if (do.segments) {
			usr <- par("usr")
			segments(x0 = v.lower[[i]],
			         x1 = v.upper[[i]],
			         y0 = h,
			         y1 = h,
			         lty = 3)
			segments(x0 = c(v.lower[[i]], v.upper[[i]]),
			         x1 = c(v.lower[[i]], v.upper[[i]]),
			         y0 = usr[3L],
			         y1 = rep.int(h, 2L),
			         lty = 3)
			text(x = mean(usr[1:2]),
			     y = h,
			     labels = sprintf("%.3g%%", 100 * level),
			     pos = 3,
			     offset = 0.1)
		}
	}

	invisible(NULL)
}
