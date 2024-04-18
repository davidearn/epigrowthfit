##' Negate Formula Terms
##'
##' Negates formula terms, simplifying \dQuote{double negatives}.
##'
##' @param x a call, symbol, or atomic constant.
##'
##' @return
##' If \code{x} is a call to \code{-} with one argument \code{y},
##' then \code{y}.  Otherwise, a call to \code{-} with \code{x} as
##' its only argument.
##'
##' @examples
##' x <- quote(x)
##' minus.x <- call("-", x)
##' stopifnot(identical(negate(x), minus.x),
##'           identical(negate(minus.x), x))

negate <-
function(x) {
	if (is.call(x) && identical(x[[1L]], quote(`-`)) && length(x) == 2L)
		x[[2L]]
	else call("-", x)
}

##' Split and Unsplit Formula Terms
##'
##' Split calls to \code{+} and \code{-} using \code{split_terms}.
##' Perform the inverse operation using \code{unsplit_terms}.
##'
##' @param x a call, symbol, or atomic constant.
##' @param l a list of calls, symbols, and atomic constants.
##'
##' @details
##' Calls to \code{(} are replaced with their first argument,
##' hence \code{x} and \code{unsplit_terms(split_terms(x))} need
##' not be identical.
##'
##' If \code{x} is a formula, then \code{split_terms(x)} operates
##' on its right hand side and \code{unsplit_terms(split_terms(x))}
##' reproduces \code{x[[length(x)]]}, not \code{x}.
##'
##' @return
##' For \code{split_terms}:
##'
##' A list of calls, symbols, and atomic constants,
##' none of which is a call to \code{+} or \code{-}.
##'
##' For \code{unsplit_terms}:
##'
##' A call, symbol, or atomic constant.
##'
##' @examples
##' x <- quote(1 + a * b - b + (c | d) + (0 + e | f))
##' l <- split_terms(x)
##' y <- unsplit_terms(l)
##' stopifnot(identical(x, y))

split_terms <-
function(x) {
	if (!(is.call(x) || is.symbol(x) || (is.atomic(x) && length(x) == 1L)))
		stop(gettextf("'%s' is not a call, symbol, or atomic constant", "x"),
		     domain = NA)
	if (inherits(x, "formula"))
		x <- x[[length(x)]]
	if (!is.call(x))
		list(x)
	else if (identical(x[[1L]], quote(`(`)))
		Recall(x[[2L]])
	else if (identical(x[[1L]], quote(`+`))) {
		if (length(x) == 2L)
			Recall(x[[2L]])
		else c(Recall(x[[2L]]), Recall(x[[3L]]))
	}
	else if (identical(x[[1L]], quote(`-`))) {
		if (length(x) == 2L)
			lapply(Recall(x[[2L]]), negate)
		else c(Recall(x[[2L]]), lapply(Recall(x[[3L]]), negate))
	}
	else list(x)
}

unsplit_terms <-
function(l) {
	if (!is.list(l))
		stop(gettextf("'%s' is not a list", "l"),
		     domain = NA)
	if (length(l) == 0L)
		return(NULL)
	x <- l[[1L]]
	while (is.call(x) && identical(x[[1L]], quote(`(`)))
		x <- x[[2L]]
	for (y in l[-1L]) {
		while (is.call(y) && identical(y[[1L]], quote(`(`)))
			y <- y[[2L]]
		x <-
			if (is.call(y) && (identical(y[[1L]], quote(`+`)) || identical(y[[1L]], quote(`-`))) && length(y) == 2L)
				as.call(list(y[[1L]], x, y[[2L]]))
			else call("+", x, y)
	}
	x
}

##' Split Fixed and Random Effect Terms
##'
##' Constructs from a mixed effects model formula the corresponding
##' fixed effects model formula and a list of all random effect terms.
##'
##' @param x a formula.
##'
##' @return
##' A list with elements:
##' \item{fixed}{a formula.}
##' \item{random}{a list of calls to \code{|}.}
##'
##' @examples
##' split_effects(y ~ 0 + x + (1 | f) + (a | g))

split_effects <-
function(x) {
	if (!inherits(x, "formula"))
		stop(gettextf("'%s' is not a formula", "x"),
		     domain = NA)
	l <- split_terms(x)
	b <- vapply(l, function(x) is.call(x) && identical(x[[1L]], quote(`|`)), FALSE)
	x[[length(x)]] <- if (all(b)) 1 else unsplit_terms(l[!b])
	list(fixed = x, random = l[b])
}

##' Split an Interaction
##'
##' Split calls to \code{:}.
##'
##' @param x a call, symbol, or atomic constant.
##'
##' @return
##' A list of calls, symbols, and atomic constants,
##' none of which is a call to \code{:}.
##'
##' @examples
##' split_interaction(quote(a:b:I(f:g):sort(h)))

split_interaction <-
function(x) {
	if (!(is.call(x) || is.symbol(x) || (is.atomic(x) && length(x) == 1L)))
		stop(gettextf("'%s' is not a call, symbol, or atomic constant", "x"),
		     domain = NA)
	if (is.call(x) && identical(x[[1L]], quote(`:`)))
		c(Recall(x[[2L]]), Recall(x[[3L]]))
	else list(x)
}

##' Replace | with + in Formula Terms
##'
##' Replaces the first element of \dQuote{first order} calls to \code{|}
##' with the symbol \code{+}, in the right hand side of a formula.
##'
##' @param x a formula.
##'
##' @examples
##' gsub_bar_plus(~x + (1 | f))

gsub_bar_plus <-
function(x) {
	if (!inherits(x, "formula"))
		stop(gettextf("'%s' is not a formula", "x"),
		     domain = NA)
	l <- split_effects(x)
	x <- l[[ "fixed"]]
	l <- l[["random"]]
	if (length(l) == 0L)
		return(x)
	l <- lapply(l, `[[<-`, 1L, quote(`+`))
	x[[length(x)]] <- unsplit_terms(c(list(x[[length(x)]]), l))
	x
}

##' Expand and Simplify Formula Terms
##'
##' Expands and simplifies calls to \code{+} and \code{-} using
##' \code{\link{terms.formula}(simplify = TRUE)} machinery.
##' Rules outlined under \code{\link{formula}} are extended to handle
##' \dQuote{first order} calls to \code{|}.
##'
##' @param x a call, symbol, or atomic constant.
##'
##' @details
##' If \code{x} is a formula, then the right hand side is expanded
##' and simplified.
##'
##' @return
##' A call, symbol, or atomic constant.
##'
##' @examples
##' simplify_terms(~0 + x * y - y)
##' simplify_terms(~0 + x * y - y + (1 | f/g))
##' simplify_terms(~0 + x * y - y + (1 | f/g) + (a | f) + (0 + b | f:g))

simplify_terms <-
function(x) {
	if (!(is.call(x) || is.symbol(x) || (is.atomic(x) && length(x) == 1L)))
		stop(gettextf("'%s' is not a call, symbol, or atomic constant", "x"),
		     domain = NA)
	if (inherits(x, "formula")) {
		x[[length(x)]] <- Recall(x[[length(x)]])
		return(x)
	}
	if (!is.call(x))
		return(x)
	l <- split_terms(x)
	b <- vapply(l, function(x) is.call(x) && identical(x[[1L]], quote(`|`)), FALSE)
	if (all(b))
		x <- NULL
	else {
		x <- terms(eval(call("~", unsplit_terms(l[!b]))), simplify = TRUE)[[2L]]
		if (!any(b))
			return(x)
	}
	Recall. <- sys.function(0L)
	expand <-
	function(x) {
		lhs <- Recall.(x[[2L]])
		rhs <- split_terms(Recall.(x[[3L]]))
		lapply(rhs, function(x) call("|", lhs, x))
	}
	condense <-
	function(l) {
		lhs <- lapply(l, `[[`, 2L)
		rhs <- l[[1L]][[3L]]
		call("|", Recall.(unsplit_terms(lhs)), rhs)
	}
	l <- unlist(lapply(l[b], expand), recursive = FALSE)
	g <- vapply(l, function(x) deparse(x[[3L]]), "")
	unsplit_terms(c(if (!is.null(x)) list(x),
	                tapply(l, g, condense, simplify = FALSE)))
}
