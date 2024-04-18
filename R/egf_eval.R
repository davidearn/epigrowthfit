##' Nonstandard Evaluation
##'
##' Utilities for obtaining index vectors from \code{subset}, \code{select},
##' and \code{order} arguments and character or expression vectors from
##' \code{label} arguments.
##'
##' @param expr
##'   a \link{symbol} or \link{call} to be evaluated in an environment
##'   constructed from \code{data} and \code{enclos} or, in the case of
##'   \code{egf_eval_select}, from
##'   \code{`names<-`(as.list(seq_len(ncol(data))), names(data))}
##'   and \code{enclos}.
##'   Alternatively, an \R{} object of a different type to be used in
##'   place of the hypothetical result of evaluation.
##' @param data
##'   a data frame.
##' @param enclos
##'   an environment to enclose \code{data}.
##'
##' @details
##' \code{subset} must evaluate to a valid index vector for
##' \code{seq_len(nrow(data))}.
##'
##' \code{select} must evaluate to a valid index vector for
##' \code{seq_len(ncol(data))}.
##'
##' \code{order} must evaluate to a permutation of
##' \code{seq_len(nrow(data))}.
##'
##' \code{label} must evaluate to a character vector or expression vector
##' of length 1 or length \code{nrow(data)} or otherwise to an \R{} object
##' coercible via \code{as.character} to such a vector.
##'
##' @return
##' Let \code{val} be the result of evaluating \code{expr}.
##'
##' \code{egf_eval_subset} returns an increasing integer vector indexing
##' rows of \code{data}, obtained after some processing of \code{val}.
##'
##' \code{egf_eval_select} returns
##' \code{match(names(data[val]), names(data), 0L)},
##' with zeros (if any) deleted.
##'
##' \code{egf_eval_order} returns \code{val}.
##'
##' \code{egf_eval_label} returns \code{val} or \code{as.character(val)},
##' repeated to length \code{nrow(data)}.
##'
##' @examples
##' year <- 2021L
##' data <- data.frame(month = sample(month.abb, 20L, replace = TRUE),
##'                    day = sample(30L, 20L, replace = TRUE))
##'
##' subset <- quote(grepl("^J", month) & day < 16L)
##' egf_eval_subset(subset, data)
##'
##' select <- quote(-day)
##' egf_eval_select(select, data)
##'
##' order <- quote(order(month, day))
##' egf_eval_order(order, data)
##'
##' label <- quote(sprintf("%04d-%02d-%02d", year, match(month, month.abb, 0L), day))
##' egf_eval_label(label, data)

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
