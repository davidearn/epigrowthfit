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
#' get_flag("test", "logspace_diff")
#'
#' @noRd
get_flag <- function(type = c("curve", "family", "prior", "test"), enum) {
  type <- match.arg(type)
  enum_all <- switch(type,
    curve = c("exponential", "subexponential", "gompertz", "logistic", "richards"),
    family = c("pois", "nbinom"),
    prior = c("norm", "lkj", "wishart", "invwishart"),
    test = c("list_of_vectors_t", "is_na", "is_finite", "logspace_diff", "mvlgamma", "log_diag_LLT", "dlkj", "dwishart", "dinvwishart", "dpois_robust", "rnbinom_robust", "eval_log_curve_exponential", "eval_log_curve_subexponential", "eval_log_curve_gompertz", "eval_log_curve_logistic", "eval_log_curve_richards", "logspace_add_baseline", "logspace_add_offsets", "eval_log_rt_subexponential", "eval_log_rt_gompertz", "eval_log_rt_logistic", "eval_log_rt_richards")
  )
  match(enum, enum_all, 0L) - 1L
}
