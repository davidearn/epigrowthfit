egf_model <-
function(curve = c("logistic", "richards",
                   "exponential", "subexponential", "gompertz"),
         excess = FALSE,
         family = c("nbinom", "pois"),
         day_of_week = FALSE) {
	stopifnot(is_true_or_false(excess), is_flag(day_of_week))
	day_of_week <- as.integer(day_of_week)
	day_of_week <- (day_of_week > 0L) * (1L + (day_of_week - 1L) %% 7L)
	res <- list(curve = match.arg(curve), excess = excess,
	            family = match.arg(family), day_of_week = day_of_week)
	class(res) <- "egf_model"
	res
}

egf_control <-
function(optimizer = egf_optimizer(),
         inner_optimizer = egf_inner_optimizer(),
         trace = FALSE,
         profile = FALSE,
         sparse_X = FALSE,
         omp_num_threads = getOption("egf.cores", 1L)) {
	if (inherits(inner_optimizer, "egf_inner_optimizer")) {
		inner_optimizer <- list(inner_optimizer)
	}
	stopifnot(inherits(optimizer, "egf_optimizer"),
	          is.list(inner_optimizer),
	          length(inner_optimizer) > 0L,
	          vapply(inner_optimizer, inherits, FALSE, "egf_inner_optimizer"),
	          is_flag(trace),
	          is_true_or_false(profile),
	          is_true_or_false(sparse_X),
	          is_number(omp_num_threads, "positive", integer = TRUE))
	trace <- as.integer(min(2, max(0, trace))) # coercion to '0:2'
	omp_num_threads <- as.integer(omp_num_threads)

	res <- list(optimizer = optimizer,
	            inner_optimizer = inner_optimizer,
	            trace = trace,
	            profile = profile,
	            sparse_X = sparse_X,
	            omp_num_threads = omp_num_threads)
	class(res) <- "egf_control"
	res
}

egf_optimizer <-
function(f = nlminb, args = list(), control = list()) {
	optimizer <- f
	stopifnot(is.list(args), is.list(control))
	if (identical(f, optim)) {
		if (is.null(args$method)) {
			args$method <- "BFGS"
			warning1("'optim' argument 'method' not specified, ",
			         "using ", dQuote(args$method), ".")
		} else {
			args$method <- match.arg(args$method, eval(formals(optim)$method))
		}
	} else if (identical(f, nlminb)) {
		f <-
		function(par, fn, gr, control, ...) {
			res <- nlminb(start = par, objective = fn, gradient = gr, control = control, ...)
			res["value"] <- res["objective"]
			res["objective"] <- NULL
			res
		}
	} else if (identical(f, nlm)) {
		f <-
		function(par, fn, gr, control, ...) {
			res <- nlm(f = structure(fn, gradient = gr), p = par, ...)
			m <- match(c("estimate", "minimum", "code"), names(res), 0L)
			names(res)[m] <- c("par", "value", "convergence")
			res["message"] <- list(NULL)
			res
		}
	} else {
		stopifnot(is.function(f),
		          length(nf <- names(formals(f))) >= 4L,
		          nf[1:3] != "...",
		          "control" %in% nf[-(1:3)])
		e <- quote(f(c(1, 1), function(x) sum(x^2), function(x) 2 * x))
		f_out <- tryCatch(eval(e),
		                  error = function(cond) {
		                  	stop("Unable to validate 'f' due to ",
		                  	     "following error in test ",
		                  	     sQuote(deparse1(e)), ":\n\n",
		                  	     conditionMessage(cond))
		                  })
		required <- c("par", "value", "convergence", "message")
		if (!(is.list(f_out) && all(required %in% names(f_out)))) {
			stop("'f' must return a list with elements ",
			     paste(sQuote(required), collapse = ", "),
			     " but _does not_ for test ", sQuote(deparse1(e)), ".")
		}
		f <-
		function(par, fn, gr, control, ...)
			optimizer(par, fn, gr, control = control, ...)
	}
	if (!is.null(names(args))) {
		reserved <- c("par", "fn", "gr", "control", "...",
		              names(formals(optimizer))[1:3])
		args[reserved] <- NULL
	}

	res <- list(f = f, args = args, control = control)
	class(res) <- "egf_optimizer"
	res
}

egf_inner_optimizer <-
function(f = newton, args = list(), control = list()) {
	stopifnot(is.list(args), is.list(control))
	if (identical(f, newton)) {
		method <- "newton"
		if (!is.null(names(args))) {
			reserved <- c("par", "fn", "gr", "he", "env", "...")
			args <- args[setdiff(names(args), reserved)]
		}
		if (!"trace" %in% names(args)) {
			args[["trace"]] <- 0L
		}
		control <- args
	} else if (identical(f, optim)) {
		if (is.null(args$method)) {
			method <- "BFGS"
			warning("'optim' argument 'method' not specified, ",
			        "using ", dQuote(method), ".")
		} else {
			method <- match.arg(args$method, eval(formals(optim)$method))
		}
	} else {
		stop("'f' is currently restricted to 'TMB::newton' and 'stats::optim'.")
	}

	res <- list(method = method, control = control)
	class(res) <- "egf_inner_optimizer"
	res
}

egf_parallel <-
function(method = c("serial", "multicore", "snow"),
         outfile = "",
         cores = getOption("egf.cores", 1L),
         args = list(),
         cl = NULL) {
	method <- match.arg(method)
	stopifnot(is_string(outfile))
	if (method == "serial") {
		cores <- 1L
		args <- list()
		cl <- NULL
	} else if (method == "multicore" || (method == "snow" && is.null(cl))) {
		stopifnot(is_number(cores, "positive", integer = TRUE), is.list(args))
		cores <- as.integer(cores)
		if (method == "multicore") {
			args$mc.cores <- cores
		} else {
			args$names <- cores
			args$outfile <- outfile
		}
	} else {
		stopifnot(inherits(cl, "SOCKcluster"))
		cores <- length(cl)
		args <- list()
	}

	res <- list(method = method, outfile = outfile,
	            cores = cores, args = args, cl = cl)
	class(res) <- "egf_parallel"
	res
}

egf_plot_control <-
function(window, data, predict, asymptote,
         box, axis, title, tdoubling, heat) {
	res <-
	list(window =
	         list(col = alpha("#DDCC77", 0.25),
	              border = NA),
	     data =
	         list(main = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD"),
	              short = list(),
	              long = list()),
	     predict =
	         list(estimate = list(col = "#44AA99", lwd = 2),
	              ci = list(border = NA)),
	     asymptote =
	         list(lty = "dotted", lwd = 2),
	     box =
	         list(bty = "l"),
	     axis =
	         list(x = list(gap.axis = 0), y = list()),
	     title =
	         list(main = list(adj = 0, xpd = NA),
	              sub = list(mgp = c(0, 0, 0), xpd = NA),
	              xlab = list(xpd = NA),
	              ylab = list(xpd = NA),
	              plab = list(col.lab = "white")),
	     tdoubling =
	         list(legend = list(),
	              estimate = list(),
	              ci = list()),
	     heat =
	         list(pal =
	                  list(colors = c("#364B9A", "#4A7BB7", "#6EA6CD",
	                                  "#98CAE1", "#C2E4EF", "#EAECCC",
	                                  "#FEDA8B", "#FDB366", "#F67E4B",
	                                  "#DD3D2D", "#A50026"),
	                       bias = 1,
	                       space = "rgb",
	                       interpolate = "linear"),
	              bg =
	                  list(col = "black",
	                       border = NA),
	              ul =
	                  list(col = alpha("black", 0.5),
	                       border = NA)))

	## Multi-assign from 'value' into 'default',
	## recursively up to one level of nesting
	rmerge <-
	function(from, into, recursive = FALSE) {
		if (is.null(from)) {
			if (recursive) {
				into[] <- list(NULL)
				return(into)
			}
			return(NULL)
		}
		if (!is.list(from) || length(from) == 0L || is.null(names(from))) {
			return(into)
		}
		nms <- unique(names(from))
		if (recursive) {
			into[nms] <- Map(rmerge,
			                 from = from[nms],
			                 into = into[nms],
			                 recursive = FALSE)
		} else {
			into[nms] <- from[nms]
		}
		into
	}

	recursive <- c("data", "predict", "axis", "title", "tdoubling", "heat")
	nms <- names(match.call()[-1L])
	if (length(nms)) {
		res[nms] <- Map(rmerge,
		                from = mget(nms, mode = "list",
		                            ifnotfound = list(list()),
		                            inherits = FALSE),
		                into = res[nms],
		                recursive = nms %in% recursive)
	}

	## Some default values are conditional on supplied values
	## FIXME: ugly ...
	for (s in c("short", "long")) {
		if (is.list(res$data[[s]])) {
			nms <- setdiff(names(res$data$main), names(res$data[[s]]))
			res$data[[s]][nms] <- res$data$main[nms]
		}
	}
	if (is.list(res$predict$ci) && is.null(res$predict$ci$col)) {
		res$predict$ci$col <- alpha(res$predict$estimate$col, 0.4)
	}
	adj <- res$title$main$adj
	if (is.numeric(adj) && length(adj) == 1L && is.finite(adj)) {
		if (is.list(res$title$sub) && is.null(res$title$sub$adj)) {
			res$title$sub$adj <- adj
		}
		if (is.list(res$tdoubling$legend) && is.null(res$tdoubling$legend$adj)) {
			res$tdoubling$legend$adj <- if (adj > 0.5) 0 else 1
		}
	}

	class(res) <- "egf_plot_control"
	res
}
