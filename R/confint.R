confint.egf <-
function(object,
         parm,
         level = 0.95,
         top = egf_top(object),
         link = TRUE,
         method = c("wald", "profile", "uniroot"),
         parallel = egf_parallel(),
         trace = FALSE,
         grid_len = 12,
         interval_scale = 7,
         subset = NULL,
         select = NULL,
         ...) {
	stopifnot(is_number_in_interval(level, 0, 1, "()"), is_true_or_false(link))
	method <- match.arg(method)
	elu <- c("estimate", "lower", "upper")

	top. <- egf_top(object)
	top <- unique(match.arg(top, top., several.ok = TRUE))

	frame_windows <- model.frame(object, "windows")
	frame_combined <- model.frame(object, "combined")
	subset <- egf_eval_subset(subset, frame_combined, parent.frame())
	select <- egf_eval_select(select, frame_combined, baseenv())

	if (method == "wald") {
		fo <- fitted(object, top = top, se = TRUE,
		             subset = subset, select = select)
		res <- confint(fo, level = level, link = link)

	}
	else if (method == "profile") {
		po <- profile(object,
		              level = level + min(0.01, 0.1 * (1 - level)),
		              top = top,
		              parallel = parallel,
		              trace = trace,
		              grid_len = grid_len,
		              subset = subset,
		              select = select)
		res <- confint(po, level = level, link = link)
		res[["linear_combination"]] <- NULL
		attr(res, "A") <- attr(res, "x") <- NULL

	}
	else { # "uniroot"
		stopifnot(inherits(parallel, "egf_parallel"),
		          is_true_or_false(trace),
		          is_number(interval_scale, "positive"))
		n <- sum(!object$random)

		l <- egf_preprofile(object, subset = subset, top = top)
		Y <- l$Y
		A <- l$A

		m <- nrow(A)
		a <- lapply(seq_len(m), function(i) A[i, ])

		## y := nll_restricted - nll_minimum = 0.5 * deviance
		target <- 0.5 * qchisq(level, df = 1)
		sd.range <- interval_scale
		nomp <- object$control$omp_num_threads

		do_uniroot <-
		function(i, a) {
			if (trace)
				cat(sprintf("Computing confidence interval %d of %d...\n",
				            i, m))
			res <- TMB::tmbroot(obj, lincomb = a, target = target,
			                    sd.range = sd.range, trace = FALSE)
			names(res) <- c("lower", "upper")
			res
		}

		if (parallel$method == "snow") {
			environment(do_uniroot) <- .GlobalEnv

			## Reconstruct list of arguments to 'MakeADFun'
			## from object internals for retaping
			args <- egf_tmb_remake_args(object$tmb_out, par = object$best)

			## Retrieve path to shared object for loading
			dll <- .dll

			cl <- parallel$cl
			if (is.null(cl)) {
				cl <- do.call(makePSOCKcluster, parallel$args)
				on.exit(stopCluster(cl), add = TRUE)
			}
			vars <- c("dll", "nomp", "args", "trace", "m", "target", "sd.range")
			clusterExport(cl, varlist = vars, envir = environment())
			clusterEvalQ(cl, {
				dyn.load(dll)
				if (TMB::openmp(n = NULL) > 0L)
					TMB::openmp(n = nomp)
				obj <- do.call(TMB::MakeADFun, args)
			})
			res <- clusterMap(cl, do_uniroot, i = seq_len(m), a = a)
		}
		else {
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
			obj <- object$tmb_out
			res <- switch(parallel$method,
			              multicore = do.call(mcMap, c(list(f = do_uniroot, i = seq_len(m), a = a), parallel$args)),
			              serial = Map(do_uniroot, i = seq_len(m), a = a))
		}

		res <- data.frame(top = rep(factor(top, levels = top.),
		                            each = length(subset)),
		                  ts = frame_windows$ts[subset],
		                  window = frame_windows$window[subset],
		                  estimate = as.double(A %*% object$best[!object$random]),
		                  do.call(rbind, res),
		                  frame_combined[subset, select, drop = FALSE],
		                  row.names = NULL,
		                  check.names = FALSE,
		                  stringsAsFactors = FALSE)
		res[elu] <- res[elu] + as.double(Y)

		if (!link) {
			f <- lapply(egf_link_extract(levels(res$top)), egf_link_match,
			            inverse = TRUE)
			res[elu] <- in_place_ragged_apply(res[elu], res$top, f = f)
			levels(res$top) <- egf_link_remove(levels(res$top))
		}
	} # "uniroot"

	row.names(res) <- NULL
	attr(res, "method") <- method
	attr(res, "level") <- level
	class(res) <- c("confint.egf", oldClass(res))
	res
}

plot.confint.egf <-
function(x,
         per_plot = 12L,
         subset = NULL,
         order = NULL,
         label = NULL,
         main = NULL,
         ...) {
	stopifnot(is_number(per_plot, "positive", integer = TRUE))
	per_plot <- as.integer(per_plot)

	subset <- egf_eval_subset(subset, x, parent.frame())
	if (length(subset) == 0L)
		stop(gettextf("'%s' is empty; nothing to plot", "subset"),
		     domain = NA)
	order <- egf_eval_order(order, x, parent.frame())
	label <- egf_eval_label(label, x, parent.frame())
	subset <- order[order %in% subset]

	a <- attributes(x)
	if (is.null(main))
		main <- sprintf("%.3g%% confidence intervals by fitting window",
		                100 * a$level)
	if (is.null(label))
		label <- as.character(x[["window"]])

	nx <- c("top", "ts", "window", "estimate", "lower", "upper")
	x <- x[nx]
	x$label <- label
	x <- x[subset, , drop = FALSE]
	x$top <- factor(x$top)

	plot.confint.egf.bars(x, per_plot = per_plot, main = main)
}

plot.confint.egf.bars <-
function(x, per_plot, main) {
	x <- split(x, x$top)
	for (xlab in names(x)) { # loop over parameters
		data <- x[[xlab]]
		n <- nrow(data)
		argna <- lapply(data[c("lower", "upper")], is.na)
		xlim <- range(data[c("estimate", "lower", "upper")], na.rm = TRUE)

		i <- 0L
		while (i < n) { # loop over plots
			k <- i + seq_len(min(per_plot, n - i))
			plot.new()
			plot.window(xlim = xlim, ylim = c(per_plot + 1, 0),
			            xaxs = "r", yaxs = "i")
			gp <- par(c("usr", "cex", "mar"))
			sfcex <- get_sfcex(x$label,
			                   target = 0.95 * max(0, gp$mar[2L] - 0.25),
			                   units = "lines")
			abline(v = axTicks(side = 1), lty = 3, col = "grey75")
			segments(x0 = replace(data$lower[k], argna$lower[k], gp$usr[1L]),
			         x1 = data$estimate[k],
			         y0 = seq_along(k),
			         y1 = seq_along(k),
			         lty = 1 + as.double(argna$lower[k]),
			         lwd = 2 - as.double(argna$lower[k]))
			segments(x0 = data$estimate[k],
			         x1 = replace(data$upper[k], argna$upper[k], gp$usr[2L]),
			         y0 = seq_along(k),
			         y1 = seq_along(k),
			         lty = 1 + as.double(argna$upper[k]),
			         lwd = 2 - as.double(argna$upper[k]))
			points(x = data$estimate[k],
			       y = seq_along(k),
			       pch = 21,
			       bg = "grey80")
			box()
			axis(side = 1)
			mtext(text = data$label[k],
			      side = 2,
			      line = 0.25,
			      at = seq_along(k),
			      las = 1,
			      adj = 1,
			      padj = 0.5,
			      cex = gp$cex * sfcex)
			title(xlab = xlab)
			title(main, adj = 0)
			i <- i + per_plot
		} # loop over plots
	} # loop over parameters

	invisible(NULL)
}
