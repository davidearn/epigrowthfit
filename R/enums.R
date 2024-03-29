##' Get flags
##'
##' Returns the underlying integer value of an enumerator in the package's
##' C++ template.
##'
##' @param enum
##'   A character vectors listing names of enumerators of type \code{type}.
##' @param type
##'   A character string specifying an enumerated type.
##'
##' @details
##' The source file defining \code{egf_get_flag} is kept synchronized with
##' the package's C++ template using \R script \file{utils/update_enums.R}.
##'
##' @return
##' An integer.
##'
##' @examples
##' egf_get_flag("exponential", "curve")
##' egf_get_flag("pois", "family")
##' egf_get_flag("norm", "prior")

egf_get_flag <-
function(enum, type) {
	enum_all <- switch(type, .__ARGS__.stop(gettextf("invalid '%s'", "type"), domain = NA))
	match(enum, enum_all, 0L) - 1L
}
