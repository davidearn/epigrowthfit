#' Get top level nonlinear model parameter names
#'
#' Retrieves the names used internally for top level nonlinear model
#' parameters.
#'
#' @param object
#'   An \R object specifying a top level nonlinear model,
#'   or \code{\link{NULL}}.
#' @param link
#'   A \link{logical} flag. If \code{TRUE},
#'   then \code{"link(name)"} is returned instead of \code{"name"}.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link{character} vector giving the subset of names relevant to
#' \code{object}, or the complete set if \code{object = \link{NULL}}.
#'
#' @examples
#' egf_get_names_top(NULL, link = FALSE)
#' egf_get_names_top(NULL, link = TRUE)
#'
#' model <- egf_model()
#' egf_get_names_top(model, link = FALSE)
#' egf_get_names_top(model, link = TRUE)
#'
#' @export
egf_get_names_top <- function(object, ...) {
  UseMethod("egf_get_names_top", object)
}

#' @rdname egf_get_names_top
#' @export
egf_get_names_top.default <- function(object, link = TRUE, ...) {
  stopifnot(is.null(object))
  names_top <- c("r", "alpha", "c0", "tinfl", "K",
                 "p", "a", "b", "disp", paste0("w", 1:6))
  if (link) {
    return(egf_link_add(names_top))
  }
  names_top
}

#' @rdname egf_get_names_top
#' @export
egf_get_names_top.egf_model <- function(object, link = TRUE, ...) {
  names_top <- switch(object$curve,
    exponential    = c("r", "c0"),
    subexponential = c("alpha", "c0", "p"),
    gompertz       = c("alpha", "tinfl", "K"),
    logistic       = c("r", "tinfl", "K"),
    richards       = c("r", "tinfl", "K", "a")
  )
  if (object$excess) {
    names_top <- c(names_top, "b")
  }
  if (object$family == "nbinom") {
    names_top <- c(names_top, "disp")
  }
  if (object$day_of_week > 0L) {
    names_top <- c(names_top, paste0("w", 1:6))
  }
  if (link) {
    return(egf_link_add(names_top))
  }
  names_top
}

#' @rdname egf_get_names_top
#' @export
egf_get_names_top.egf <- function(object, link = TRUE, ...) {
  egf_get_names_top(object$model, link = link)
}

#' @export
egf_get_names_top.tmb_data <- function(object, link = TRUE, ...) {
  names_top <- levels(object$X_info$par)
  if (link) {
    return(names_top)
  }
  egf_link_remove(names_top)
}

#' Check for random effects
#'
#' Determines whether an object specifies a random effects model.
#'
#' @param object An \code{"\link{egf}"} object.
#'
#' @return
#' \code{TRUE} or \code{FALSE}.
#'
#' @examples
#' e <- new.env()
#' e$data <- list()
#' object <- list(tmb_out = list(env = e))
#' class(object) <- "egf"
#'
#' e$data$Z <- matrix(numeric(9L), 3L, 3L) # at least one column
#' egf_has_random(object)
#'
#' e$data$Z <- matrix(numeric(0L), 3L, 0L) # zero columns
#' egf_has_random(object)
#'
#' @export
egf_has_random <- function(object) {
  stopifnot(inherits(object, "egf"))
  ncol(object$tmb_out$env$data$Z) > 0L
}

#' Construct the combined model frame
#'
#' Joins in a single \link[=data.frame]{data frame} all mixed effects
#' \link[=model.frame]{model frames} and further variables originally
#' indicated by \code{\link{egf}} argument \code{append}.
#'
#' @param object An \code{"\link{egf}"} object.
#'
#' @details
#' If a variable name occurs in multiple mixed effects model frames,
#' then only the first instance is retained. Except in unusual cases
#' (possible only if model formulae have different formula environments),
#' all instances of a variable name are identical, and no information is lost.
#'
#' Since the data frames being combined each correspond rowwise
#' to \code{object$frame_windows}, so does the result.
#'
#' @return
#' A \link[=data.frame]{data frame} combining
#' (in the sense of \code{\link{cbind}})
#' all data frames in \link{list} \code{object$frame_parameters}
#' and data frame \code{object$frame_append}.
#'
#' @examples
#' object <- list(
#'   frame_parameters = list(a = data.frame(x = seq_len(9L)), b = data.frame(y = rnorm(9L))),
#'   frame_append = data.frame(x = rnorm(9L), z = 0)
#' )
#' class(object) <- "egf"
#' egf_combine_frames(object)
#'
#' @export
egf_combine_frames <- function(object) {
  stopifnot(inherits(object, "egf"))
  l <- c(unname(object$frame_parameters), list(object$frame_append))
  res <- do.call(cbind, l)
  res[!duplicated(names(res))]
}
