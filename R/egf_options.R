egf_model <-
function(curve = c("logistic", "richards",
                   "exponential", "subexponential", "gompertz"),
         excess = FALSE,
         family = c("nbinom", "pois"),
         day_of_week = FALSE) {
	stopifnot(isTrueFalse(excess), isInteger(day_of_week))
	day_of_week <- as.integer(day_of_week)
	day_of_week <- (day_of_week > 0L) * (1L + (day_of_week - 1L) %% 7L)
	ans <- list(curve = match.arg(curve), excess = excess,
	            family = match.arg(family), day_of_week = day_of_week)
	class(ans) <- "egf_model"
	ans
}

egf_control <-
function(outer_optimizer = egf_optimizer(nlminb),
         inner_optimizer = egf_optimizer(newton),
         trace = FALSE,
         profile = FALSE,
         sparse_X = FALSE,
         omp_num_threads = getOption("egf.cores", 1L)) {
	stopifnot(inherits(outer_optimizer, "egf_optimizer"),
	          inherits(inner_optimizer, "egf_optimizer"),
	          isInteger(trace),
	          isTrueFalse(profile),
	          isTrueFalse(sparse_X),
	          isInteger(omp_num_threads), omp_num_threads >= 1)
	trace <- as.integer(min(2, max(0, trace))) # coercion to 0:2
	omp_num_threads <- as.integer(omp_num_threads)

	ans <- list(outer_optimizer = outer_optimizer,
	            inner_optimizer = inner_optimizer,
	            trace = trace,
	            profile = profile,
	            sparse_X = sparse_X,
	            omp_num_threads = omp_num_threads)
	class(ans) <- "egf_control"
	ans
}

egf_optimizer <-
function(f = nlminb, args = list(), control = list()) {
	stopifnot(is.list(args), is.list(control))
	optimizer <- f
	if (identical(f, newton)) {
		f <- .f <- # avoid false positive NOTE: "multiple local ..."
		function(par, fn, gr, he = function(par) optimHess(par, fn, gr),
		         control = control, ...) {
			ans <- newton(par = par, fn = fn, gr = gr, he = he,
			              control = control, ...)
			ans[c("convergence", "message")] <- list(0L, NULL)
			ans
		}
		if (!any(names(args) == "trace"))
			args[["trace"]] <- 0L
		control <- list()
		attr(f, "method") <- "newton"
	}
	else if (identical(f, optim)) {
		args[["method"]] <-
			if (is.null(args[["method"]]))
				"BFGS"
			else match.arg(args[["method"]], eval(formals(optim)[["method"]]))
		attr(f, "method") <- "optim"
	}
	else if (identical(f, nlminb)) {
		f <-
		function(par, fn, gr, control, ...) {
			ans <- nlminb(start = par, objective = fn, gradient = gr,
			              control = control, ...)
			m <- match("objective", names(ans), 0L)
			names(ans)[m] <- "value"
			ans
		}
		attr(f, "method") <- "nlminb"
	}
	else if (identical(f, nlm)) {
		f <-
		function(par, fn, gr, control, ...) {
			ans <- nlm(f = `attr<-`(fn, "gradient", gr), p = par, ...)
			m <- match(c("estimate", "minimum", "code"), names(ans), 0L)
			names(ans)[m] <- c("par", "value", "convergence")
			ans["message"] <- list(NULL)
			ans
		}
		attr(f, "method") <- "nlm"
	}
	else {
		stopifnot(is.function(f),
		          length(nf <- names(formals(f))) >= 4L,
		          all(nf[1:3] != "..."),
		          any(nf[-(1:3)] == "control"))
		f.val <- f(c(1, 1), function(x) sum(x * x), function(x) 2 * x)
		required <- c("par", "value", "convergence", "message")
		if (!(is.list(f.val) && all(match(required, names(f.val), 0L))))
			stop(gettextf("'%s' must return a list with elements %s",
			              "f", deparse(required)),
			     domain = NA)
		f <-
		function(par, fn, gr, control, ...)
			optimizer(par, fn, gr, control = control, ...)
		attr(f, "method") <- ""
	}
	if (!is.null(names(args))) {
		reserved <- c("par", "fn", "gr", "he", "control", "...",
		              names(formals(optimizer))[1:3])
		args <- args[match(names(args), reserved, 0L) == 0L]
	}
	ans <- list(f = f, args = args, control = control)
	class(ans) <- "egf_optimizer"
	ans
}

egf_parallel <-
function(method = c("serial", "multicore", "snow"),
         outfile = "",
         cores = getOption("egf.cores", 1L),
         args = list(),
         cl = NULL) {
	method <- match.arg(method)
	stopifnot(isString(outfile))
	if (method == "serial") {
		cores <- 1L
		args <- list()
		cl <- NULL
	}
	else if (method == "multicore" || (method == "snow" && is.null(cl))) {
		stopifnot(isInteger(cores), cores >= 1, is.list(args))
		cores <- as.integer(cores)
		if (method == "multicore")
			args[["mc.cores"]] <- cores
		else {
			args[["names"]] <- cores
			args[["outfile"]] <- outfile
		}
	}
	else {
		stopifnot(inherits(cl, "SOCKcluster"))
		cores <- length(cl)
		args <- list()
	}
	ans <- list(method = method, outfile = outfile,
	            cores = cores, args = args, cl = cl)
	class(ans) <- "egf_parallel"
	ans
}

.defaults <-
	list(window =
	         list(col = "#DDCC773F", # alpha("#DDCC77", 0.25)
	              border = NA),
	     data =
	         list(main = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD"),
	              short = list(),
	              long = list()),
	     predict =
	         list(value = list(col = "#44AA99", lwd = 2),
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
	              ylab = list(xpd = NA)),
	     doubling =
	         list(legend = list(),
	              value = list(),
	              ci = list()))

egf_control_plot <-
function(window, data, predict, asymptote, box, axis, title, doubling) {
	ans <- .defaults

	## Multi-assign from 'value' into 'default', recursively
	rmerge <-
	function(from, into, depth = 0L) {
		if (is.null(from))
			return(if (depth > 0L) { into[] <- list(NULL); into })
		if (!is.list(from) || length(from) == 0L || is.null(names(from)))
			return(into)
		nms <- unique(names(from))
		into[nms] <-
			if (depth > 0L)
				.mapply(rmerge, list(from = from[nms], into = into[nms]), list(depth = depth - 1L))
			else from[nms]
		into
	}

	recursive <- c("data", "predict", "axis", "title", "doubling")
	nms <- names(match.call())[-1L]
	if (length(nms))
		ans[nms] <-
			.mapply(rmerge,
			        list(from = mget(nms, mode = "list",
			                         ifnotfound = list(list()),
			                         inherits = FALSE),
			             into = ans[nms],
			             depth = as.integer(match(nms, recursive, 0L) > 0L)),
			        NULL)

	## Some default values are conditional on supplied values
	## FIXME: ugly ...
	for (s in c("short", "long"))
		if (is.list(ans[["data"]][[s]])) {
			nms <- setdiff(names(ans[["data"]][["main"]]),
			               names(ans[["data"]][[s]]))
			ans[["data"]][[s]][nms] <- ans[["data"]][["main"]][nms]
		}
	if (is.list(ans[["predict"]][["ci"]]) &&
	    is.null(ans[["predict"]][["ci"]][["col"]]))
		ans[["predict"]][["ci"]][["col"]] <-
			alpha(ans[["predict"]][["value"]][["col"]], 0.4)
	if (!is.null(title.main.adj <- ans[["title"]][["main"]][["adj"]])) {
	if (is.list(ans[["title"]][["sub"]]) &&
	    is.null(ans[["title"]][["sub"]][["adj"]]))
		ans[["title"]][["sub"]][["adj"]] <-
			title.main.adj
	if (is.list(ans[["doubling"]][["legend"]]) &&
	    is.null(ans[["doubling"]][["legend"]][["adj"]]))
		ans[["doubling"]][["legend"]][["adj"]] <-
			if (title.main.adj > 0.5) 0 else 1
	}
	class(ans) <- "egf_control_plot"
	ans
}
