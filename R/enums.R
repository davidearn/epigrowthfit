#' Get flags
#'
#' Returns the underlying integer value of an enumerator in the package's
#' C++ template.
#'
#' @param type
#'   A \link{character} string. The name of an enumerated type.
#' @param enum
#'   A \link{character} string. The name of an enumerator of type \code{type}.
#'
#' @details
#' The source file defining \code{get_flag} is kept synchronized with
#' the package's C++ template using R script \file{utils/update_enums.R}.
#'
#' @return
#' An \link{integer}.
#'
#' @examples
#' get_flag("curve", "exponential")
#' get_flag("family", "pois")
#' get_flag("prior", "norm")
#'
#' @noRd
get_flag <- function(type, enum) {
  curve_names <- c("exponential", "subexponential", "gompertz", "logistic", "richards")
  family_names <- c("pois", "nbinom")
  prior_names <- c("norm")
  switch(type,
    curve = match(enum, curve_names) - 1L,
    family = match(enum, family_names) - 1L,
    prior = match(enum, prior_names) - 1L,
    stop("Type not implemented.")
  )
}
