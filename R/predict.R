#' Compute predicted values
#'
#' Computes predicted values of log cumulative incidence, log
#' interval incidence, and log instantaneous exponential growth
#' rate given user-specified time points.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param time
#'   A numeric vector listing increasing time points in days since
#'   the start of the fitting window specified by `subset`.
#' @param varname
#'   A character vector listing one or more variables for which
#'   predicted values are desired.
#' @param subset
#'   A named list of atomic vectors of length 1 specifying exactly
#'   one level for each factor in `object$frame` (and thus a unique
#'   fitting window). Use the default (`NULL`) if and only if
#'   `object$frame` has no factors.
#' @param se
#'   A logical scalar. If `TRUE`, then standard errors on predicted
#'   values are also reported.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Elements of `subset` (if non-`NULL`) must have the form
#' `factor_name = level_name`, where `factor_name` is
#' the name of a factor in `object$frame` and `level_name`
#' is an element of `levels(object$frame$factor_name)`.
#'
#' @return
#' A list inheriting from class `"egf_predict"`,
#' containing data frames (one for each `varname`)
#' with numeric variables `time`, `estimate`, and
#' `se` (if `se = TRUE`).
#'
#' If `"log_int_inc" %in% varname`, then, for `i > 1`,
#' `log_int_inc$estimate[i]` stores log incidence
#' during the interval from `time[i-1]` to `time[i]`.
#'
#' @seealso [confint.egf_predict()]
#' @export
#' @importFrom TMB MakeADFun sdreport
predict.egf <- function(object,
                        time,
                        varname = c("log_int_inc", "log_cum_inc", "log_rt"),
                        subset = NULL,
                        se = FALSE,
                        ...) {
  varname <- unique(match.arg(varname, several.ok = TRUE))
  index_levels <- levels(object$index)
  i1 <- which(!duplicated(object$index))
  refdate <- object$frame[i1, 1L]
  if (length(object$frame) > 2L || !is.null(subset)) {
    frame_red <- object$frame[i1, -(1:2), drop = FALSE]
    stop_if_not(
      is.list(subset),
      length(subset) > 0L,
      !is.null(names(subset)),
      m = "`subset` must be a named list or NULL."
    )
    stop_if_not(
      length(subset) == length(frame_red),
      identical(sort(names(subset)), sort(names(frame_red))),
      vapply(subset, is.atomic, logical(1L)),
      lengths(subset) == 1L,
      mapply(`%in%`, subset, lapply(frame_red[names(subset)], levels)),
      m = paste0(
        "`subset` must specify exactly one level\n",
        "for each factor in `object$frame`."
      )
    )
    subset <- lapply(subset, as.character)
    w <- Reduce(`&`, Map(`==`, frame_red[names(subset)], subset))
    stop_if_not(
      any(w),
      m = "`subset` does not match any fitting windows."
    )
    index_levels <- index_levels[w]
    i1 <- i1[w]
    refdate <- refdate[w]
  }
  stop_if_not(
    is.numeric(time),
    length(time) >= 2L,
    m = "`time` must be numeric and have length 2 or greater."
  )
  stop_if_not(
    !anyNA(time),
    m = "`time` must not have missing values."
  )
  stop_if_not(
    diff(time) > 0,
    m = "`time` must be increasing."
  )
  stop_if_not_tf(se)
  tr <- range(object$tmb_args$data$t[object$index %in% index_levels])
  warn_if_not(
    time >= tr[1L],
    time <= tr[2L],
    m = paste0(
      "There are elements of `time` outside of\n",
      "the fitting window specified by `subset`."
    )
  )

  ## Create the data objects needed to run prediction code
  ## in C++ template
  sparse_X_flag <- Xs <- Xd <- Z <- NULL # R CMD check
  object$tmb_args$data <- within(object$tmb_args$data, {
    predict_flag <- 1L
    predict_lci_lii_lrt_flag <- 1L * (c("log_int_inc", "log_cum_inc", "log_rt") %in% varname)
    se_flag <- 1L * se
    t_new <- time
    if (sparse_X_flag == 1L) {
      ## sparseMatrix() doesn't recycle like matrix()
      Xs_new <- do.call(rbind, rep.int(list(Xs[i1[1L], , drop = FALSE]), length(t_new)))
      Xd_new <- Xd
    } else {
      Xd_new <- matrix(Xd[i1[1L], , drop = TRUE], nrow = length(t_new), ncol = ncol(Xd), byrow = TRUE)
      Xs_new <- Xs
    }
    if (has_random(object)) {
      Z_new <- do.call(rbind, rep.int(list(Z[i1[1L], , drop = FALSE]), length(t_new)))
    } else {
      Z_new <- Z
    }
  })

  tmb_out_new <- do.call(MakeADFun, object$tmb_args)
  varname_new <- sprintf("%s_new", varname)
  if (se) {
    tmb_out_new$fn(object$par[object$nonrandom])
    out <- split_sdreport(sdreport(tmb_out_new))[varname_new]
  } else {
    out <- lapply(tmb_out_new$report(object$par)[varname_new],
                  function(x) data.frame(estimate = x))
  }
  if ("log_int_inc" %in% varname) {
    out$log_int_inc_new <- rbind(NA_real_, out$log_int_inc_new)
  }
  out <- lapply(out, function(d) cbind(time, d))
  names(out) <- sub("_new$", "", names(out))
  attr(out, "subset") <- subset
  attr(out, "refdate") <- refdate
  class(out) <- c("egf_predict", "list")
  out
}

#' Confidence bands on predicted values
#'
#' Computes confidence bands on predicted log cumulative incidence,
#' predicted log interval incidence, and predicted log instantaneous
#' exponential growth rate, assuming asymptotic normality.
#'
#' @param object
#'   An `"egf_predict"` object returned by [predict.egf()].
#'   Must supply standard errors on log predicted incidence.
#' @param log
#'   A logical scalar. If `FALSE`, then log scale predicted
#'   values are inverse log transformed.
#' @inheritParams confint.egf_profile
#'
#' @details
#' Confidence bands on log scale predicted values are computed as
#' `estimate + c(-1, 1) * q * se`, where `q = qnorm(0.5 * (1 + level))`,
#' where `estimate` and `se` are the log predicted values and the
#' approximate (delta method) standard errors.
#'
#' @return
#' If `log = TRUE`, then `object` but with variable `se` in each listed
#' data frame replaced with two variables `lower` and `upper` supplying
#' confidence bands on log predicted values.
#'
#' Otherwise, the same list but with variables `estimate`, `lower`,
#' and `upper` inverse log transformed in each listed data frame.
#' In this case, the prefix `"log_"` is stripped from the list names.
#'
#' @export
#' @importFrom stats qnorm
confint.egf_predict <- function(object, parm, level = 0.95,
                                log = TRUE, ...) {
  stop_if_not(
    "se" %in% names(object[[1L]]),
    m = paste0(
      "`object` must supply standard errors.\n",
      "Repeat `predict()` with argument `se = TRUE`."
    )
  )
  stop_if_not_in_0_1(level)

  q <- qnorm(0.5 * (1 + level))
  f <- if (log) identity else exp
  g <- function(d) {
    elu <- d$estimate + outer(d$se, c(0, -1, 1) * q)
    colnames(elu) <- c("estimate", "lower", "upper")
    data.frame(time = d$time, f(elu))
  }
  out <- lapply(object, g)
  if (!log) {
    names(out) <- sub("^log_", "", names(out))
  }
  attr(out, "subset") <- attr(object, "subset")
  attr(out, "refdate") <- attr(object, "refdate")
  attr(out, "level") <- level
  out
}
