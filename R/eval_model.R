#' \loadmathjax
#' Evaluate expected cumulative incidence
#'
#' @description
#' Evaluates expected cumulative incidence at desired time points
#' conditional on a vector of model parameters.
#'
#' @param time A numeric vector listing time points in days since
#'   a reference date.
#' @param curve One of `"exponential"`, `"logistic"`, and `"richards"`,
#'   indicating a model of expected cumulative incidence.
#' @param include_baseline A logical scalar. If `TRUE`, then
#'   the model of expected cumulative incidence will include
#'   a linear baseline.
#' @param theta A named numeric vector listing positive values
#'   for all relevant model parameters:
#'
#'   \describe{
#'     \item{`r`}{Initial exponential growth rate, expressed per day.}
#'     \item{`c0`}{Expected cumulative incidence on the reference date.
#'       Used only if `curve = "exponential"`.
#'     }
#'     \item{`K`}{Expected epidemic final size.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`thalf`}{Time at which the epidemic is expected to attain
#'       half its final size, expressed as a number of days since the
#'       reference date.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`p`}{Richards shape parameter.
#'       Used only if `curve = "richards"`.
#'     }
#'     \item{`b`}{Baseline linear growth rate, expressed per day.
#'       Used only if `include_baseline = TRUE`.
#'     }
#'   }
#'
#' @return
#' A numeric vector with length `length(time)` whose `i`th
#' element is expected cumulative incidence at `time[i]`.
#'
#' @details
#' A full description of the models of expected cumulative
#' incidence specified by `curve` and `include_baseline`
#' can be found in the package vignette, accessible with
#' `vignette("epigrowthfit-vignette")`.
#'
#' @export
eval_model <- function(time, curve, include_baseline = FALSE, theta) {
  if (!is.numeric(time)) {
    stop("`time` must be numeric.")
  } else if (isTRUE(any(diff(time) < 0))) {
    stop("`time[!is.na(time)]` must be increasing.")
  }
  if (!is.character(curve) || length(curve) != 1 ||
      !curve %in% c("exponential", "logistic", "richards")) {
    stop("`curve` must be an element of ",
         "`c(\"exponential\", \"logistic\", \"richards\")`.")
  }
  if (!is.logical(include_baseline) || length(include_baseline) != 1 ||
      is.na(include_baseline)) {
    stop("`include_baseline` must be `TRUE` or `FALSE`.")
  }
  if (!is.numeric(theta) || is.null(names(theta))) {
    stop("`theta` must be a named numeric vector.")
  }
  par <- switch(curve,
    exponential = c("r", "c0"),
    logistic    = c("r", "K", "thalf"),
    richards    = c("r", "K", "thalf", "p")
  )
  if (include_baseline) {
    par <- c(par, "b")
  }
  if (!all(par %in% names(theta))) {
    stop("`theta` must specify all model parameters.")
  } else if (!all(is.finite(theta[par])) || !all(theta[par] > 0)) {
    stop("`theta` must specify positive values for all model parameters.")
  }

  with(as.list(theta[par]), {
    x <- switch(curve,
      exponential = c0 * exp(r * time),
      logistic = K / (1 + exp(-r * (time - thalf))),
      richards = K / (1 + (2^p - 1) * exp(-r * p * (time - thalf)))^(1 / p)
    )
    if (include_baseline) {
      x <- x + b * time
    }
    x
  })
}
