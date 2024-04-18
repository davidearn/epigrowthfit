##' Get Flags
##'
##' Returns the underlying integer value of an enumerator in the package's
##' \proglang{C++} template.
##'
##' @param enum
##'   a character vector listing names of enumerators of type \code{type}.
##' @param type
##'   a character string specifying an enumerated type.
##'
##' @details
##' The source file defining \code{egf_get_flag} is kept synchronized with
##' the package's \proglang{C++} template using the package \file{Makevars}.
##'
##' @return
##' An integer vector.
##'
##' @examples
##' egf_get_flag("exponential", "curve")
##' egf_get_flag("pois", "family")
##' egf_get_flag("norm", "prior")

egf_get_flag <-
function(enum, type) {
	enum.all <- switch(type, .__ARGS__.stop(gettextf("invalid '%s'", "type"), domain = NA))
	match(enum, enum.all, 0L) - 1L
}
