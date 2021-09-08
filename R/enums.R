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
#' The source file defining \code{egf_get_flag} is kept synchronized with
#' the package's C++ template using R script \file{utils/update_enums.R}.
#'
#' @return
#' An \link{integer}.
#'
#' @examples
#' egf_get_flag("curve", "exponential")
#' egf_get_flag("family", "pois")
#' egf_get_flag("prior", "norm")
#'
#' @noRd
egf_get_flag <- function(type = c("curve", "family", "prior", "test"), enum) {
  type <- match.arg(type)
  enum_all <- switch(type,
    curve = c("exponential", "subexponential", "gompertz", "logistic", "richards"),
    family = c("pois", "nbinom"),
    prior = c("norm", "lkj", "wishart", "invwishart")
  )
  match(enum, enum_all, 0L) - 1L
}
