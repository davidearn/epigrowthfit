#' \loadmathjax
#' Evaluate expected cumulative incidence
#'
#' @description
#' Evaluates expected cumulative incidence at desired time points
#' conditional on a vector of model parameters.
#'
#' @param time A numeric vector listing time points in days.
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
#'     \item{`c0`}{Expected cumulative incidence at `time = 0`.
#'       Used only if `curve = "exponential"`.
#'     }
#'     \item{`K`}{Expected epidemic final size.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`thalf`}{Time in days at which the epidemic is expected
#'       to attain half its final size.
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
  check(time,
    what = "numeric",
    "`time` must be numeric."
  )
  check(curve,
    what = "character",
    len = 1,
    opt = c("exponential", "logistic", "richards"),
    "`curve` must be one of ",
    "\"exponential\", \"logistic\", \"richards\"."
  )
  check(include_baseline,
    what = "logical",
    len = 1,
    opt = c(TRUE, FALSE),
    "`include_baseline` must be TRUE or FALSE."
  )
  par <- switch(curve,
    exponential = c("r", "c0"),
    logistic    = c("r", "K", "thalf"),
    richards    = c("r", "K", "thalf", "p")
  )
  if (include_baseline) {
    par <- c(par, "b")
  }
  check(theta,
    what = "numeric",
    no = function(x) is.null(names(x)),
    "`theta` must be a named numeric vector."
  )
  check(theta[par],
    yes = function(x) all(is.finite(x) & x > 0),
    "`theta` must specify finite positive values for:\n",
    paste(par, collapse = ", ")
  )

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
