egf_eval_subset <-
function(expr, data, enclos = parent.frame()) {
	stopifnot(is.data.frame(data))
	n <- nrow(data)
	sln <- seq_len(n)
	if (is.null(expr))
		return(sln)
	names(sln) <- row.names(data)
	subset <-
		if (is.symbol(expr) || is.call(expr))
			eval(expr, data, enclos)
		else expr
	ans <- sln[subset]
	names(ans) <- NULL
	sort(unique(ans[!is.na(ans)]))
}

egf_eval_select <-
function(expr, data, enclos = baseenv()) {
	stopifnot(is.data.frame(data))
	if (is.null(expr))
		return(integer(0L))
	nms <- names(data)
	select <-
		if (is.symbol(expr) || is.call(expr))
			eval(expr, `names<-`(as.list(seq_along(data)), nms), enclos)
		else expr
	ans <- match(names(data[select]), nms, 0L)
	ans[ans > 0L]
}

egf_eval_order <-
function(expr, data, enclos = parent.frame()) {
	stopifnot(is.data.frame(data))
	n <- nrow(data)
	sln <- seq_len(n)
	if (is.null(expr))
		return(sln)
	order <-
		if (is.symbol(expr) || is.call(expr))
			eval(expr, data, enclos)
		else expr
	stopifnot(is.integer(order), length(order) == n, all(sort(order) == sln))
	order
}

egf_eval_label <-
function(expr, data, enclos = parent.frame()) {
	stopifnot(is.data.frame(data))
	if (is.null(expr))
		return(NULL)
	n <- nrow(data)
	label <-
		if (is.symbol(expr) || is.call(expr))
			eval(expr, data, enclos)
		else expr
	if (!(is.character(label) || is.expression(label)))
		label <- as.character(label)
	stopifnot(any(length(label) == c(1L, n)))
	rep_len(label, n)
}
