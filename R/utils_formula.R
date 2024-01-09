##' Negate formula terms
##'
##' Negates formula terms, taking some care to simplify double negatives.
##'
##' @param x A \link{call}, \link{name}, or \link{atomic} scalar.
##'
##' @return
##' If \code{x} is a \link{call} to the \code{`-`} operator with
##' one argument \code{y}, then \code{y}. Otherwise, a call to
##' the \code{`-`} operator with \code{x} as its only argument.
##'
##' @examples
##' x <- quote(x)
##' minus_x <- call("-", x)
##' stopifnot(identical(negate(x), minus_x))
##' stopifnot(identical(negate(minus_x), x))

negate <-
function(x) {
	if (is.call(x) && x[[1L]] == "-" && length(x) == 2L) {
		return(x[[2L]])
	}
	call("-", x)
}

##' Split and unsplit formula terms
##'
##' Split nested \link{call}s to binary operators \code{`+`}
##' and \code{`-`} with \code{split_terms}.
##' Perform the inverse operation with \code{unsplit_terms}.
##'
##' @param x
##'   A \link{call}, \link{name}, or \link{atomic} scalar.
##' @param l
##'   A \link{list} of \link{call}s, \link{name}s, and \link{atomic} scalars.
##'
##' @details
##' Semantically extraneous calls to \code{`(`},
##' as in the expression \code{(x) + (y)} are stripped when splitting.
##' Calls to \code{`|`} are wrapped in \code{`(`} when unsplitting.
##' For these reasons, it is possible for
##' \code{\link{identical}(x, unsplit_terms(split_terms(x)))}
##' to return \code{FALSE}.
##'
##' If \code{x} is a formula, then the right hand side is split.
##' In this case, \code{unsplit_terms(split_terms(x))} reproduces
##' \code{x[[\link{length}(x)]]}, not \code{x}.
##'
##' @return
##' A \link{list} of \link{call}s, \link{name}s, and \link{atomic} scalars.
##'
##' @examples
##' x <- quote(1 + a * b - b + (c | d) + (0 + e | f))
##' l <- split_terms(x)
##' y <- unsplit_terms(l)
##' stopifnot(identical(x, y))

split_terms <-
function(x) {
	if (inherits(x, "formula")) {
		x <- x[[length(x)]]
	}
	if (is.name(x) || (is.atomic(x) && length(x) == 1L)) {
		return(list(x))
	}
	stopifnot(is.call(x))
	if (x[[1L]] == "(") {
		return(split_terms(x[[2L]]))
	}
	if (x[[1L]] == "+") {
		if (length(x) == 2L) {
			return(split_terms(x[[2L]]))
		} else {
			return(c(split_terms(x[[2L]]), split_terms(x[[3L]])))
		}
	}
	if (x[[1L]] == "-") {
		if (length(x) == 2L) {
			return(lapply(split_terms(x[[2L]]), negate))
		} else {
			return(c(split_terms(x[[2L]]), split_terms(x[-2L])))
		}
	}
	list(x)
}

unsplit_terms <-
function(l) {
	stopifnot(is.list(l))
	if (length(l) == 0L) {
		return(NULL)
	}
	is_pm <-
	function(x)
		is.call(x) && (x[[1L]] == "+" || x[[1L]] == "-") && length(x) == 2L
	is_bar <-
	function(x)
		is.call(x) && x[[1L]] == "|"
	x <- l[[1L]]
	if (is_bar(x)) {
		x <- call("(", x)
	}
	for (i in seq_along(l)[-1L]) {
		y <- l[[i]]
		if (is_pm(y)) {
			x <- as.call(list(y[[1L]], x, y[[2L]]))
		} else {
			if (is_bar(y)) {
				y <- call("(", y)
			}
			x <- call("+", x, y)
		}
	}
	x
}

##' Split fixed and random effects terms
##'
##' Retrieves from a mixed effects model \link{formula} the corresponding
##' fixed effects model formula and a \link{list} of all random effects terms.
##'
##' @param x A \link{formula}.
##'
##' @return
##' A \link{list} with elements:
##' \item{fixed}{
##'   A fixed effects model \link{formula}.
##' }
##' \item{random}{
##'   A \link{list} of random effects terms,
##'   i.e., \link{call}s to the \code{`|`} operator.
##' }
##'
##' @examples
##' split_effects(y ~ 0 + x + (1 | f) + (a | g))

split_effects <-
function(x) {
	stopifnot(inherits(x, "formula"))
	l <- split_terms(x)
	is_bar <-
	function(x)
		is.call(x) && x[[1L]] == "|"
	l_is_bar <- vapply(l, is_bar, FALSE)
	x[[length(x)]] <- if (all(l_is_bar)) 1 else unsplit_terms(l[!l_is_bar])
	list(fixed = x, random = l[l_is_bar])
}

##' Split an interaction
##'
##' Recursively constructs a \link{list} of arguments to binary operator
##' \code{`:`} from a nested \link{call} to \code{`:`}, excluding arguments
##' that are themselves \link{call}s to \code{`:`}.
##'
##' @param x A \link{call}, \link{name}, or \link{atomic} scalar.
##'
##' @return
##' A \link{list} of \link{call}s, \link{name}s, and \link{atomic} scalars.
##'
##' @examples
##' x <- quote(a:b:I(f:g):sort(h))
##' split_interaction(x)

split_interaction <-
function(x) {
	if (is.name(x) || (is.atomic(x) && length(x) == 1L)) {
		return(list(x))
	}
	if (is.call(x)) {
		if (x[[1L]] == ":") {
			return(do.call(c, lapply(x[-1L], split_interaction)))
		} else {
			return(list(x))
		}
	}
	stop("'x' must be a call, name, or atomic scalar.")
}

##' Replace `|` with `+` in formula terms
##'
##' Replaces the \code{`|`} operator in mixed effects model \link{formula}e
##' with the \code{`+`} operator, enabling construction of model frames using
##' \code{\link{model.frame}} machinery.
##'
##' @param x A \link{formula}.
##'
##' @return
##' \code{x} with \code{`|`} in terms of the form \code{(tt | g)} replaced
##' with \code{`+`}.
##'
##' @examples
##' gsub_bar_plus(~x + (1 | f))

gsub_bar_plus <-
function(x) {
	stopifnot(inherits(x, "formula"))
	l <- split_effects(x)
	if (length(l$random) == 0L) {
		return(x)
	}
	m <- length(x)
	l$random <- lapply(l$random, `[[<-`, 1L, as.name("+"))
	x[[m]] <- unsplit_terms(c(l$fixed[[m]], l$random))
	x
}

##' Expand and simplify formula terms
##'
##' Expands and simplifies nested \link{call}s to binary operators
##' \code{`+`} and \code{`-`} using \code{\link{terms}(simplify = TRUE)}
##' machinery. Rules outlined under \code{\link{formula}} are extended
##' to address handling of calls to binary operator \code{`|`}.
##'
##' @param x
##'   A \link{call}, \link{name}, or \link{atomic} scalar.
##'
##' @details
##' \code{x} is split into a \link{list} of constituent terms.
##' Terms that _are not_ \link{call}s to \code{`|`} are regrouped,
##' and the resulting expression is expanded and simplified
##' using \code{\link{terms.formula}(simplify = TRUE)}.
##' Terms that _are_ \link{call}s to \code{`|`} are regrouped,
##' and the resulting expression is expanded and simplified
##' following extended rules (see below).
##' Finally, the two expressions are merged, yielding the final result.
##'
##' For example, consider the expression \code{w + (x * y | f) + (z | f/g)}.
##' The first component of the final result is simply \code{w}.
##' The second component is obtained from the subexpression
##' \code{(x * y | f) + (z | f/g)}.
##' First, the arguments of each \link{call} to \code{`|`} are
##' expanded and simplified using \code{\link{terms.formula}(simplify = TRUE)},
##' producing \code{(x + y + x:y | f) + (z | f + f:g)}.
##' Second, left hand expressions are distributed to right hand terms,
##' producing \code{(x + y + x:y | f) + (z | f) + (z | f:g)}.
##' Third, left hand expressions with matching right hand terms are merged,
##' producing \code{(x + y + x:y + z | f) + (z | f:g)}.
##' Finally, left hand expressions are simplified,
##' again using \code{\link{terms.formula}(simplify = TRUE)},
##' producing \code{(x + y + z + x:y | f) + (z | f:g)}.
##' (In this case, the effect is merely a permutation of terms
##' according to their order.) Hence the final result is
##' \code{w + (x + y + z + x:y | f) + (z | f:g)}.
##'
##' If \code{x} is a \link{formula}, then the left and right hand
##' expressions are expanded and simplified separately.
##'
##' @return
##' The \code{call}, \link{name}, or \link{atomic} scalar resulting from
##' expansion and simplification of the constituent terms of \code{x}.
##'
##' @examples
##' simplify_terms(~0 + x * y - y)
##' simplify_terms(~0 + x * y - y + (1 | f/g))
##' simplify_terms(~0 + x * y - y + (1 | f/g) + (a | f) + (0 + b | f:g))

simplify_terms <-
function(x) {
	if (inherits(x, "formula")) {
		x[-1L] <- lapply(x[-1L], simplify_terms)
		return(x)
	}
	if (is.name(x) || (is.atomic(x) && length(x) == 1L)) {
		return(x)
	}
	stopifnot(is.call(x))
	l <- split_terms(x)
	is_bar <-
	function(x)
		is.call(x) && x[[1L]] == "|"
	l_is_bar <- vapply(l, is_bar, FALSE)
	if (all(l_is_bar)) {
		no_bar <- NULL
	} else {
		no_bar <- terms(as.formula(call("~", unsplit_terms(l[!l_is_bar]))), simplify = TRUE)[[2L]]
		if (!any(l_is_bar)) {
			return(no_bar)
		}
	}
	expand_bar <-
	function(x) {
		lhs <- simplify_terms(x[[2L]])
		rhs <- split_terms(simplify_terms(x[[3L]]))
		lapply(rhs, function(x) call("|", lhs, x))
	}
	bar <- do.call(c, lapply(l[l_is_bar], expand_bar))
	rhs_deparsed <- vapply(bar, function(x) deparse1(x[[3L]]), "")
	merge_bars <-
	function(l) {
		lhs <- lapply(l, `[[`, 2L)
		rhs <- l[[1L]][[3L]]
		call("|", simplify_terms(unsplit_terms(lhs)), rhs)
	}
	bar <- tapply(bar, rhs_deparsed, merge_bars, simplify = FALSE)
	unsplit_terms(c(no_bar, bar))
}
