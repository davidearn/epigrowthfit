confint.egf <-
function(object, parm, level = 0.95,
         A = seq_along(par),
         method = c("wald", "profile", "uniroot"),
         scale = 7,
         parallel = egf_parallel(),
         trace = FALSE,
         top = egf_top(object),
         subset = NULL,
         select = NULL,
         class = FALSE,
         link = TRUE,
         random = FALSE,
         ...) {
	stopifnot(isNumber(level), level > 0, level < 1,
	          isTrueFalse(class), isTrueFalse(link), isTrueFalse(random))
	method <- match.arg(method)

	par <- coef(object)
	which <- NULL

	if (m.A <- is.null(A)) {
		top. <- egf_top(object)
		top <- unique(match.arg(top, top., several.ok = TRUE))

		frame <- model.frame(object, "combined")
		subset <- egf_eval_subset(subset, frame, parent.frame())
		select <- egf_eval_select(select, frame, baseenv())
		frame <- frame[subset, select, drop = FALSE]

		info <- model.frame(object, "windows")[subset, c("ts", "window")]
		ts     <- info[["ts"    ]]
		window <- info[["window"]]

		A <- egf_make_A(object, top = top, subset = subset)
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

	h <- 0.5 * (1 - level)
	p <- c(h, 1 - h)
	q <- qnorm(p)
	percent <- formatp(p, 3L)
	if (class)
		class <- m.A

	if (method == "wald" && m.A && random) {
		rpt <- egf_adreport(object)
		i <- names(rpt[["value"]]) == "Y"
		value <- rpt[["value"]][i]
		se    <- rpt[["sd"   ]][i]
		dim(value) <- dim(se) <- rpt[["env"]][["ADreportDims"]][["Y"]]
		dimnames(value) <- dimnames(se) <- list(NULL, top.)
		value <- as.vector(value[subset, top, drop = FALSE])
		se    <- as.vector(se   [subset, top, drop = FALSE])
		ans <- wald(value, se, level)
	}
	else if (method == "wald") {
		rpt <- egf_adreport(object)
		P. <- rpt[["par.fixed"]]
		V. <- rpt[["cov.fixed"]]
		if (!is.null(which)) {
			value <- P.[which]
			se    <- sqrt(diag(V., names = FALSE)[which])
		}
		else {
			A. <- A[, -1L, drop = FALSE]
			V. <- A. %*% tcrossprod(V., A.)
			value <- as.vector(A %*% c(1, P.))
			se    <- sqrt(diag(V., names = FALSE))
		}
		ans <- wald(value, se, level)
	}
	else if (method == "profile") {
		po <- profile(object, level = level, A = A,
		              parallel = parallel, trace = trace, ...)
		ans <- confint(po, level = level)
	}
	else if (method == "uniroot") {
		stopifnot(inherits(parallel, "egf_parallel"),
		          isTrueFalse(trace),
		          isNumber(scale), scale > 0)

		a <- lapply(seq_len(nrow(A)), function(i) A[i, ])

		## y := nll_restricted - nll_minimum = 0.5 * (change in deviance)
		target <- 0.5 * qchisq(level, df = 1)
		sd.range <- scale

		nomp <- object[["control"]][["omp_num_threads"]]
		nprf <- length(a)

		doRoot <-
		function(i, a) {
			if (trace)
				cat(gettextf("computing interval %d of %d ...", i, nprf),
				    "\n", sep = "")
			a[1L] + TMB::tmbroot(obj, lincomb = a[-1L], target = target,
			                     sd.range = sd.range, trace = FALSE)
		}

		if (parallel[["method"]] == "snow") {
			environment(doRoot) <- globalenv()

			## Reconstruct list of arguments to 'MakeADFun'
			## from object internals for retaping
			args <- egf_tmb_remake_args(object[["tmb_out"]], par = object[["best"]])

			## Retrieve path to shared object for loading
			dll <- .dll

			cl <- parallel[["cl"]]
			if (is.null(cl)) {
				cl <- do.call(makePSOCKcluster, parallel[["args"]])
				on.exit(stopCluster(cl), add = TRUE)
			}
			vars <- c("dll", "nomp", "nprf", "args", "trace",
			          "target", "sd.range")
			clusterExport(cl, varlist = vars, envir = environment())
			clusterEvalQ(cl, {
				dyn.load(dll)
				if (TMB::openmp(n = NULL) > 0L)
					TMB::openmp(n = nomp)
				obj <- do.call(TMB::MakeADFun, args)
			})
			ans <- clusterMap(cl, doRoot, i = seq_len(nprf), a = a)
		}
		else {
			if (given_outfile <- nzchar(parallel[["outfile"]])) {
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
			obj <- object[["tmb_out"]]
			ans <- switch(parallel[["method"]],
			              multicore = do.call(mcMap, c(list(f = doRoot, i = seq_len(nprf), a = a), parallel[["args"]])),
			              serial = Map(doRoot, i = seq_len(nprf), a = a))
		}
		ans <- matrix(unlist1(ans), length(ans), 2L, byrow = TRUE)
	}

	if (class || !m.A)
		dimnames(ans) <- list(if (!class) dimnames(A)[[1L]], percent)
	else {
		dim(ans) <- c(length(subset), length(top), 2L)
		dimnames(ans) <- list(paste(ts, window, sep = ", "), top, percent)
	}

	if (class) {
		ns <- nrow(frame)
		nt <- length(top)
		i. <- rep.int(seq_len(ns), nt)
		j. <- rep(seq_len(nt), each = ns)
		if (!link) {
			pos <- seq.int(from = 0L, by = ns, length.out = nt + 1L)
			for (j in seq_len(length(pos) - 1L)) {
				f <- egf_link_match(egf_link_extract(top[j]), inverse = TRUE)
				k <- seq.int(pos[j] + 1L, pos[j + 1L])
				ans[k, ] <- f(ans[k, ])
			}
			top <- egf_link_remove(top)
			top. <- egf_link_remove(top.)
		}
		tmp <- data.frame(top = factor(top, levels = top.)[j.],
		                  ts = ts[i.],
		                  window = window[i.],
		                  value =
		                      if (method == "wald" && m.A && random)
		                          value
		                      else as.vector(A %*% c(1, par)),
		                  ci = NA,
		                  frame[i., ],
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

plot.confint.egf <-
function(x, by = 12L,
         subset = NULL, order = NULL, label = NULL, main = NULL, ...) {
	stopifnot(isInteger(by), by >= 1)
	by <- as.integer(by)

	subset <- egf_eval_subset(subset, x, parent.frame())
	order <- egf_eval_order(order, x, parent.frame())
	subset <- order[match(order, subset, 0L) > 0L]

	x <- x[subset, , drop = FALSE]
	nrx <- nrow(x)
	ncx <- length(x)

	value <- x[["value"]]
	lower <- x[["ci"]][, 1L]
	upper <- x[["ci"]][, 2L]
	label <-
		if (is.null(label))
			as.character(x[["window"]])
		else egf_eval_label(label, x, parent.frame())

	ylim <- c(by + 1, 0)

	grp <- split(seq_len(nrx), factor(x[["top"]]))
	nms <- names(grp)
	for (g in seq_along(grp)) { # loop over parameters
		l <- grp[[g]]
		value. <- value[l]
		lower. <- lower[l]
		upper. <- upper[l]
		label. <- label[l]

		xlim <- c(min(lower., na.rm = TRUE), max(upper., na.rm = TRUE))
		xlab <- nms[g]

		pos <- c(seq.int(from = 0L, to = length(l) - 1L, by = by), length(l))
		for (j in seq_len(length(pos) - 1L)) { # loop over plots
			k <- seq.int(pos[j] + 1L, pos[j + 1L]);
			value.. <- value.[k]
			lower.. <- lower.[k]
			upper.. <- upper.[k]
			label.. <- label.[k]

			h <- as.double(seq_along(k))

			dev.hold()
			plot.new()
			plot.window(xlim = xlim, ylim = ylim, xaxs = "r", yaxs = "i")
			usr <- par("usr")
			mar <- par("mar")
			cex <- par("cex") *
				get_sfcex(label, 0.95 * max(0, mar[2L] - 0.25), "lines")
			abline(v = axTicks(side = 1L), lty = 3L, col = "grey75")
			argna <- is.na(lower..)
			segments(x0 = replace(lower.., argna, usr[1L]),
			         x1 = value..,
			         y0 = h,
			         y1 = h,
			         lty = 1L + argna,
			         lwd = 2L - argna)
			argna <- is.na(upper..)
			segments(x0 = value..,
			         x1 = replace(upper.., argna, usr[2L]),
			         y0 = h,
			         y1 = h,
			         lty = 1L + argna,
			         lwd = 2L - argna)
			points(x = value..,
			       y = h,
			       pch = 21L,
			       bg = "grey80")
			box()
			axis(side = 1L)
			mtext(text = label..,
			      side = 2L,
			      line = 0.25,
			      at = h,
			      las = 1,
			      adj = 1,
			      padj = 0.5,
			      cex = cex)
			title(main = main)
			title(xlab = xlab)
			dev.flush()
		} # loop over plots
	} # loop over parameters

	invisible(NULL)
}
