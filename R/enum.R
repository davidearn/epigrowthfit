#' Get flags
#'
#' Returns the underlying integer value of an enumerator in
#' the C++ template.
#'
#' @param type
#'   A character string. The name of an enumerated type.
#' @param enum
#'   A character string. The name of an enumerator of type `type`.
#'
#' @details
#' `get_flag()` is kept synchronized with the C++ template
#' using the R script `utils/update_enum.R`.
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
