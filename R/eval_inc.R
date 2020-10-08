#' \loadmathjax
#' Evaluate expected incidence
#'
#' @description
#' Evaluates expected cumulative and interval incidence at desired time
#' points, conditional on a phenomenological model of cumulative incidence.
#'
#' @param time A numeric vector listing (increasing) time points in days
#'   since some initial date.
#' @param curve One of `"exponential"`, `"logistic"`, and `"richards"`,
#'   indicating a phenomenological model for cumulative incidence.
#' @param include_baseline A logical scalar. If `TRUE`, then the
#'   cumulative incidence model will include a linear term \mjseqn{b t}.
#' @param theta A named numeric vector listing positive values for all
#'   relevant model parameters:
#'
#'   \describe{
#'     \item{`r`}{\mjseqn{\lbrace\,r\,\rbrace}
#'       Initial (exponential) growth rate expressed per day.
#'     }
#'     \item{`x0`}{\mjseqn{\lbrace\,x_0\,\rbrace}
#'       Expected cumulative incidence on the initial date.
#'       This is the expectation of the number of cases observed
#'       up to the initial date. Used only if `curve = "exponential"`.
#'     }
#'     \item{`K`}{\mjseqn{\lbrace\,K\,\rbrace}
#'       Expected epidemic final size. This is the expectation of the
#'       number of cases observed over the full course of the epidemic.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`thalf`}{\mjseqn{\lbrace\,t_\text{half}\,\rbrace}
#'       Expected time at which the epidemic attains half its
#'       final size, expressed as a (possibly non-integer) number
#'       of days since the initial date.
#'       Used only if `curve %in% c("logistic", "richards")`.
#'     }
#'     \item{`p`}{\mjseqn{\lbrace\,p\,\rbrace}
#'       Richards shape parameter. Used only if `curve = "richards"`.
#'     }
#'     \item{`b`}{\mjseqn{\lbrace\,b\,\rbrace}
#'       Baseline (linear) growth rate expressed per day.
#'       Used only if `include_baseline = TRUE`.
#'     }
#'   }
#'
#' @return
#' A list with elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`cum_inc`}{A numeric vector with length `length(time)`
#'     such that `cum_inc[i]` is the expected number of cases observed
#'     up to `time[i]`.
#'   }
#'   \item{`int_inc`}{A numeric vector with length `length(time)-1`,
#'     such that `int_inc[i]` is the expected number of cases observed
#'     between `time[i]` and `time[i+1]`. Equal to `diff(cum_inc)`.
#'   }
#' }
#'
#' @details
#' Equations specifying how expected cumulative incidence is computed
#' for a given model are presented in Details 1 and 2 of [egf_init()].
#'
#' @export
eval_inc <- function(time, curve, include_baseline = FALSE, theta) {
  if (missing(time)) {
    stop("Missing argument `time`.")
  } else if (!is.numeric(time)) {
    stop("`time` must be numeric.")
  } else if (isTRUE(any(diff(time) < 0))) {
    stop("`time[!is.na(time)]` must be increasing.")
  }
  if (missing(curve)) {
    stop("Missing argument `curve`.")
  } else if (!is.character(curve) || length(curve) != 1 ||
               !curve %in% c("exponential", "logistic", "richards")) {
    stop("`curve` must be an element of ",
         "`c(\"exponential\", \"logistic\", \"richards\")`.")
  }
  if (!is.logical(include_baseline) || length(include_baseline) != 1 ||
        is.na(include_baseline)) {
    stop("`include_baseline` must be `TRUE` or `FALSE`.")
  }
  if (missing(theta)) {
    stop("Missing argument `theta`.")
  } else if (!is.numeric(theta) || is.null(names(theta))) {
    stop("`theta` must be a named numeric vector.")
  }
  par <- switch(curve,
    exponential = c("r", "x0"),
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

  cum_inc <- with(as.list(theta[par]), {
    x <- switch(curve,
      exponential = x0 * exp(r * time),
      logistic = K / (1 + exp(-r * (time - thalf))),
      richards = K / (1 + (2^p - 1) * exp(-r * p * (time - thalf)))^(1 / p)
    )
    if (include_baseline) {
      x <- x + b * time
    }
    x
  })
  list(time = time, cum_inc = cum_inc, int_inc = diff(cum_inc))
}
