##' Manipulate top level nonlinear model parameter names
##'
##' Utilities for modifying strings \code{s} and
##' \code{fs = \link{sprintf}("\%s(\%s)", f, s)},
##' where \code{s} is the name used internally
##' for a top level nonlinear model parameter
##' and \code{f} is the name of the link function
##' (either \code{"log"} or \code{"logit"})
##' used internally for that parameter.
##'
##' @param s,fs Character vectors.
##'
##' @return
##' For strings \code{s}, \code{f}, and \code{fs} as described above,
##' \code{egf_link_get(s)} returns \code{f},
##' \code{egf_link_add(s)} returns \code{fs},
##' \code{egf_link_remove(fs)} returns \code{s}, and
##' \code{egf_link_extract(fs)} returns \code{f}.
##'
##' \code{NA_character_} is returned for invalid values of the argument.
##'
##' @examples
##' egf_link_get("r")
##' egf_link_add("r")
##' egf_link_remove("log(r)")
##' egf_link_extract("log(r)")
##'
##' egf_link_extract("invalid string", "r" , "log(r)")

egf_link_get <-
function(s) {
	ok <- s %in% egf_get_names_top(NULL, link = FALSE)
	s[ok] <- replace(rep.int("log", sum(ok)), s[ok] == "p", "logit")
	s[!ok] <- NA
	s
}

egf_link_add <-
function(s) {
	ok <- s %in% egf_get_names_top(NULL, link = FALSE)
	s[ok] <- sprintf("%s(%s)", egf_link_get(s[ok]), s[ok])
	s[!ok] <- NA
	s
}

egf_link_remove <-
function(fs) {
	ok <- fs %in% egf_get_names_top(NULL, link = TRUE)
	fs[ok] <- sub("^(log|logit)\\((.*)\\)$", "\\2", fs[ok])
	fs[!ok] <- NA
	fs
}

egf_link_extract <-
function(fs) {
	ok <- fs %in% egf_get_names_top(NULL, link = TRUE)
	fs[ok] <- sub("^(log|logit)\\((.*)\\)$", "\\1", fs[ok])
	fs[!ok] <- NA
	fs
}

##' Get link and inverse link functions
##'
##' Retrieve the link function named by a string, or its inverse.
##'
##' @param f
##'   A \link{character} string naming a link function.
##' @param inverse
##'   A \link{logical} flag. If \code{TRUE}, then the inverse is returned.
##'
##' @return
##' A function.
##'
##' @examples
##' egf_link_match("log")
##' egf_link_match("log", inverse = TRUE)

egf_link_match <-
function(f, inverse = FALSE) {
	if (inverse) {
		switch(f,
		       identity = identity,
		       log = exp,
		       logit = function(q) plogis(q),
		       stop("Link not implemented."))
	} else {
		switch(f,
		       identity = identity,
		       log = log,
		       logit = function(p) qlogis(p),
		       stop("Link not implemented."))
	}
}
