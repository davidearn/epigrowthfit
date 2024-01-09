egf_eval_subset <-
function(expr, data, enclos = parent.frame()) {
	stopifnot(is.data.frame(data))
	n <- nrow(data)
	if (is.null(expr)) {
		return(seq_len(n))
	}
	if (is.language(expr)) {
		subset <- eval(expr, data, enclos)
	} else {
		subset <- expr
	}
	index <- seq_len(n)
	names(index) <- row.names(data)
	subset <- index[subset]
	sort(unique(subset[!is.na(subset)]))
}

egf_eval_order <-
function(expr, data, enclos = parent.frame()) {
	stopifnot(is.data.frame(data))
	n <- nrow(data)
	if (is.null(expr)) {
		return(seq_len(n))
	}
	if (is.language(expr)) {
		order <- eval(expr, data, enclos)
	} else {
		order <- expr
	}
	stopifnot(is.numeric(order),
	          length(order) == n,
	          sort(order) == seq_len(n))
	as.integer(order)
}

egf_eval_append <-
function(expr, data, enclos = baseenv()) {
	stopifnot(is.data.frame(data))
	if (is.null(expr)) {
		return(integer(0L))
	}
	if (is.language(expr)) {
		l <- as.list(seq_along(data))
		names(l) <- names(data)
		append <- eval(expr, l, enclos)
	} else {
		append <- expr
	}
	append <- match(names(data[append]), names(data), 0L)
	append[append > 0L]
}

egf_eval_label <-
function(expr, data, enclos = parent.frame()) {
	stopifnot(is.data.frame(data))
	if (is.null(expr)) {
		return(NULL)
	}
	if (is.language(expr)) {
		label <- eval(expr, data, enclos)
	} else {
		label <- expr
	}
	label <- as.character(label)
	stopifnot(length(label) %in% c(1L, n <- nrow(data)))
	rep_len(label, n)
}
