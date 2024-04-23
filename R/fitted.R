fitted.egf <-
function(object,
         top = egf_top(object),
         subset = NULL,
         select = NULL,
         class = FALSE,
         se = FALSE,
         ...) {
	stopifnot(isTrueFalse(class), isTrueFalse(se))
	if (se)
		se <- class

	top. <- egf_top(object)
	top <- unique(match.arg(top, top., several.ok = TRUE))

	frame <- model.frame(object, "combined")
	subset <- egf_eval_subset(subset, frame, parent.frame())
	select <- egf_eval_select(select, frame, baseenv())
	frame <- frame[subset, select, drop = FALSE]

	info <- model.frame(object, "windows")[subset, c("ts", "window")]
	ts     <- info[["ts"    ]]
	window <- info[["window"]]

	if (se. <- se) {
		rpt <- egf_adreport(object)
		i <- names(rpt[["value"]]) == "Y"
		value <- rpt[["value"]][i]
		se    <- rpt[["sd"   ]][i]
		dim(value) <- dim(se) <- rpt[["env"]][["ADreportDims"]][["Y"]]
		dimnames(value) <- dimnames(se) <- list(NULL, top.)
		value <- as.vector(value[subset, top, drop = FALSE])
		se    <- as.vector(se   [subset, top, drop = FALSE])
	}
	else {
		value <- egf_report(object)[["Y"]]
		dimnames(value) <- list(NULL, top.)
		value <- as.vector(value[subset, top, drop = FALSE])
	}

	ns <- length(subset)
	nt <- length(top)

	if (class) {
		i. <- rep.int(seq_len(ns), nt)
		j. <- rep(seq_len(nt), each = ns)
		ans <- data.frame(top = factor(top, levels = top.)[j.],
		                  ts = ts[i.],
		                  window = window[i.],
		                  value = value,
		                  se = NA,
		                  frame[i., ],
		                  row.names = NULL,
		                  check.names = FALSE,
		                  stringsAsFactors = FALSE)
		ans[["se"]] <- if (se.) se
		attr(ans, "se") <- se.
		attr(ans, "ns") <- ns
		attr(ans, "nt") <- nt
		class(ans) <- c("fitted.egf", oldClass(ans))
	}
	else {
		dim(value) <- c(ns, nt)
		dimnames(value) <- list(paste(ts, window, sep = ", "), top)
		ans <- value
	}

	ans
}

fitted.egf_no_fit <-
function(object, ...) {
	call <- match.call(expand.dots = FALSE)
	call[[1L]] <- quote(fitted.egf)
	object[["best"]] <- object[["init"]]
	eval(call)
}

confint.fitted.egf <-
function(object, parm = seq_len(nrow(object)), level = 0.95,
         class = FALSE, ...) {
	stopifnot(isNumber(level), level > 0, level < 1,
	          isTrueFalse(class), attr(object, "se"))
	which <- seq_len(nrow(object))[parm]
	if (anyNA(which))
		stop(gettextf("invalid '%s'", "parm"), domain = NA)
	o. <- object[which, ]
	attr(o., "se") <- attr(o., "ns") <- attr(o., "nt") <- NULL
	class(o.) <- oldClass(o.)[oldClass(o.) != "fitted.egf"]

	h <- 0.5 * (1 - level)
	p <- c(h, 1 - h)
	q <- qnorm(p)
	percent <- formatp(p, 3L)

	ans <- wald(o.[["value"]], o.[["se"]], level)

	if (class) {
		dimnames(ans) <- list(NULL, percent)
		row.names(o.) <- NULL
		names(o.)[match("se", names(o.), 0L)] <- "ci"
		o.[["ci"]] <- ans
		attr(o., "level") <- level
		class(o.) <- c("confint.egf", oldClass(o.))
		ans <- o.
	}
	else if (!missing(parm))
		dimnames(ans) <- list(paste(o.[["top"]], o.[["ts"]], o.[["window"]], sep = ", "),
		                      percent)
	else {
		ns <- attr(object, "ns")
		nt <- attr(object, "nt")
		dim(ans) <- c(ns, nt, 2L)
		dimnames(ans) <- list(paste(o.[["ts"]][seq_len(ns)], o.[["window"]][seq_len(ns)], sep = ", "),
		                      as.character(o.[["top"]][seq.int(1L, by = ns, length.out = nt)]),
		                      percent)
	}

	ans
}
