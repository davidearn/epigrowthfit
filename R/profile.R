profile.egf <-
function(fitted,
         level = 0.95,
         which = NULL,
         A = NULL,
         top = egf_get_names_top(fitted, link = TRUE),
         parallel = egf_parallel(),
         trace = FALSE,
         grid_len = 12,
         subset = NULL,
         append = NULL,
         .subset = NULL,
         .append = NULL,
         ...) {
	stopifnot(is_number_in_interval(level, 0, 1, "()"),
	          inherits(parallel, "egf_parallel"),
	          is_true_or_false(trace),
	          is_number_in_interval(grid_len, 1, Inf, "[)"))
	n <- sum(!fitted$random)

	## If profiling user-specified elements of 'c(beta, theta)'
	if (!is.null(which)) {
		method <- "which"
		eval(bquote(stopifnot(is.numeric(which),
		                      which %in% seq_len(.(n)))))
		which <- unique(which)
		A <- sparseMatrix(i = seq_len(m),
		                  j = which,
		                  x = 1,
		                  dims = c(length(which), n))

	## If profiling user-specified linear combinations
	## of elements of 'c(beta, theta)'
	} else if (!is.null(A)) {
		method <- "A"
		if (!is(A, "dMatrix")) {
			stopifnot(is.numeric(A))
			if (is.null(dim(A))) {
				dim(A) <- c(1L, length(A))
			} else {
				stopifnot(is.matrix(A))
			}
		}
		eval(bquote(stopifnot(nrow(A) > 0L,
		                      ncol(A) == .(n),
		                      is.finite(A),
		                      rowSums(abs(A)) > 0)))

	## If profiling population fitted values
	## of top level nonlinear model parameters
	} else if (!is.null(top)) {
		method <- "top"

		names_top <- egf_get_names_top(fitted, link = TRUE)
		top <- unique(match.arg(top, names_top, several.ok = TRUE))

		frame_windows <- model.frame(fitted, "windows")
		frame_combined <- model.frame(fitted, "combined")
		subset <- if (is.null(.subset)) substitute(subset) else .subset
		subset <- egf_eval_subset(subset, frame_combined, parent.frame())
		append <- if (is.null(.append)) substitute(append) else .append
		append <- egf_eval_append(append, frame_combined, baseenv())

		l <- egf_preprofile(fitted, subset = subset, top = top)
		Y <- l$Y
		A <- l$A

	## Otherwise
	} else {
		stop("One of 'which', 'A', and 'top' must be non-NULL.")
	}

	## Covariance matrix of 'c(beta, theta)'
	V <- unclass(vcov(fitted))
	m <- nrow(A)
	if (method == "which") {
		a <- which
		h <- 0.25 * sqrt(diag(V)[which])
		s <- "name"
	} else {
		## Covariance matrix of 'A %*% c(beta, theta)'
		V <- A %*% tcrossprod(V, A)
		a <- lapply(seq_len(m), function(i) A[i, ])
		h <- 0.25 * sqrt(diag(V))
		s <- "lincomb"
	}

	## y := nll_restricted - nll_minimum = 0.5 * deviance
	ytol <- 0.5 * qchisq(level, df = 1)
	ystep <- ytol / grid_len
	nomp <- fitted$control$omp_num_threads

	do_profile <-
	function(i, a, h) {
		if (trace) {
			cat(sprintf("Computing likelihood profile %d of %d...\n", i, m))
		}
		args <- list(obj = obj, h = h, ytol = ytol, ystep = ystep)
		args[[s]] <- a
		res <- do.call(TMB::tmbprofile, args)
		i_min <- which.min(res[[2L]])
		## deviance = 2 * (nll_restricted - nll_minimum)
		res[[2L]] <- 2 * (res[[2L]] - res[i_min, 2L])
		names(res) <- c("value", "deviance")
		res[-i_min, , drop = FALSE] # 'tmbprofile' duplicates this row
	}

	if (parallel$method == "snow") {
		environment(do_profile) <- .GlobalEnv

		## Reconstruct list of arguments to 'MakeADFun' from object internals
		## for retaping
		args <- egf_tmb_remake_args(fitted$tmb_out, par = fitted$best)

		## Retrieve path to shared object for loading
		dll <- system.file("libs", TMB::dynlib("epigrowthfit"),
		                   package = "epigrowthfit", mustWork = TRUE)

		cl <- parallel$cl
		if (is.null(cl)) {
			cl <- do.call(makePSOCKcluster, parallel$args)
			on.exit(stopCluster(cl), add = TRUE)
		}
		vars <- c("dll", "nomp", "args", "trace", "m", "s", "ytol", "ystep")
		clusterExport(cl, varlist = vars, envir = environment())
		clusterEvalQ(cl, {
			dyn.load(dll)
			if (TMB::openmp(n = NULL) > 0L) {
				TMB::openmp(n = nomp)
			}
			obj <- do.call(TMB::MakeADFun, args)
		})
		res <- clusterMap(cl, do_profile, i = seq_len(m), a = a, h = h)
	} else {
		if (given_outfile <- nzchar(parallel$outfile)) {
			outfile <- file(parallel$outfile, open = "wt")
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
		obj <- fitted$tmb_out
		res <- switch(parallel$method,
		              multicore = do.call(mcMap, c(list(f = do_profile, i = seq_len(m), a = a, h = h), parallel$args)),
		              serial = Map(do_profile, i = seq_len(m), a = a, h = h))
	}

	nr <- vapply(res, nrow, 0L)
	res <- data.frame(linear_combination = rep.int(gl(m, 1L), nr),
	                  do.call(rbind, res),
	                  row.names = NULL)
	if (method == "top") {
		i <- rep.int(rep.int(subset, length(top)), nr)
		res <- data.frame(top = rep.int(rep(factor(top, levels = names_top),
		                                    each = length(subset)), nr),
		                  frame_windows[i, c("ts", "window"), drop = FALSE],
		                  res,
		                  frame_combined[i, append, drop = FALSE],
		                  row.names = NULL,
		                  check.names = FALSE,
		                  stringsAsFactors = FALSE)
		res$value <- res$value + rep.int(as.double(Y), nr)
	}
	attr(res, "A") <- A
	attr(res, "x") <- fitted$best[!fitted$random]
	attr(res, "level") <- level
	attr(res, "method") <- method
	class(res) <- c("egf_profile", oldClass(res))
	res
}

confint.egf_profile <-
function(object, parm, level = attr(object, "level"), link = TRUE, ...) {
	stopifnot(is_true_or_false(link),
	          is_number_in_interval(level, 0, attr(object, "level"), "(]"))
	q <- qchisq(level, df = 1)
	method <- attr(object, "method")

	do_solve <-
	function(d) {
		i_min <- which.min(d$deviance)
		i_left <- seq_len(i_min)
		i_right <- seq.int(i_min, nrow(d))
		estimate <- d$value[i_min]
		lower <- approx(x = d$deviance[i_left],
		                y = d$value[i_left],
		                xout = q)$y
		upper <- approx(x = d$deviance[i_right],
		                y = d$value[i_right],
		                xout = q)$y
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

	res <- do.call(rbind, by(object, object$linear_combination, do_solve,
	                         simplify = FALSE))
	if (method == "top" && !link) {
		elu <- c("estimate", "lower", "upper")
		f <- lapply(egf_link_extract(levels(res$top)), egf_link_match,
		            inverse = TRUE)
		res[elu] <- in_place_ragged_apply(res[elu], res$top, f = f)
		levels(res$top) <- egf_link_remove(levels(res$top))
	}
	res$linear_combination <- as.integer(as.character(res$linear_combination))
	row.names(res) <- NULL
	attr(res, "A") <- attr(object, "A")
	attr(res, "x") <- attr(object, "x")
	attr(res, "level") <- level
	res
}

plot.egf_profile <-
function(x, level = attr(x, "level"), sqrt = FALSE, subset = NULL, ...) {
	subset <- egf_eval_subset(substitute(subset), x, parent.frame())
	subset <- match(levels(factor(x$linear_combination[subset])),
	                levels(x$linear_combination))

	stopifnot(is_true_or_false(sqrt))
	f <- if (sqrt) base::sqrt else identity
	do_segments <-
		!sqrt &&
		is.numeric(level) &&
		length(level) > 0L &&
		!all(is.na(level))
	if (do_segments) {
		level <- level[!is.na(level)]
		stopifnot(level > 0, level <= attr(x, "level"))
		## Line segments at heights 'h' in all plots
		h <- qchisq(level, df = 1)
		## Line segment 'j' to start at 'v_lower[[i]][j]'
		## and end at 'v_upper[[i]][j]' in plot 'i'
		ci <- lapply(level, function(p) confint(x, level = p))
		v_lower <- lapply(subset, function(i) vapply(ci, `[`, 0, i, "lower"))
		v_upper <- lapply(subset, function(i) vapply(ci, `[`, 0, i, "upper"))
	}

	dots <- list(...)
	method <- attr(x, "method")
	if (is.null(dots[["main"]])) {
		if (method == "top") {
			main <- c(tapply(as.character(x$window),
			                 x$linear_combination, `[[`, 1L))[subset]
		} else {
			main <- character(length(subset))
		}
	} else {
		main <- rep_len(dots$main, length(subset))
	}
	if (is.null(dots[["xlab"]])) {
		if (method == "top") {
			xlab <- c(tapply(as.character(x$top),
			                 x$linear_combination, `[[`, 1L))[subset]
		} else {
			xlab <- paste0("linear combination ",
			               levels(x$linear_combination)[subset])
		}
	} else {
		xlab <- rep_len(dots$xlab, length(subset))
	}
	dots$main <- dots$xlab <- NULL
	if (is.null(dots[["ylab"]])) {
		dots$ylab <- "deviance"
		if (sqrt) {
			dots$ylab <- bquote(sqrt(.(dots$ylab)))
		}
	}
	if (is.null(dots[["ylim"]])) {
		dots$ylim <- c(0, f(max(x$deviance, na.rm = TRUE)))
	}

	for (i in seq_along(subset)) {
		args <- list(formula = f(deviance) ~ value,
		             data = x,
		             subset = (as.integer(x$linear_combination) == subset[i]),
		             main = main[i],
		             xlab = xlab[i])
		do.call(plot, c(args, dots))
		if (do_segments) {
			usr <- par("usr")
			segments(x0 = v_lower[[i]],
			         x1 = v_upper[[i]],
			         y0 = h,
			         y1 = h,
			         lty = 3)
			segments(x0 = c(v_lower[[i]], v_upper[[i]]),
			         x1 = c(v_lower[[i]], v_upper[[i]]),
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
