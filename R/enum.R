#' Get flags
#'
#' Returns the underlying integer value of an enumerator in the
#' C++ template.
#'
#' @param type
#'   A character string giving the name of an enumerated type.
#' @param enum
#'   A character string giving the name of an enumerator of type `type`.
#'
#' @details
#' This function is kept synchronized with the C++ template using script
#' `utils/update_flags.R`.
#'
#' @return
#' An integer.
#'
#' @examples
#' get_flag("curve", "logistic")
#' get_flag("distr", "nbinom")
#'
#' @noRd
get_flag <- function(type, enum) {
  curve_names <- c("exponential", "subexponential", "gompertz", "logistic", "richards") # GREP_FLAG
  distr_names <- c("pois", "nbinom") # GREP_FLAG
  switch(type,
    curve = match(enum, curve_names) - 1L,
    distr = match(enum, distr_names) - 1L,
    NA_integer_
  )
}
