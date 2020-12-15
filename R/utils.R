stop_if_not <- function(..., m) {
  n <- ...length()
  for (i in seq_len(n)) {
    e <- ...elt(i)
    if (!(is.logical(e) && !anyNA(e) && all(e))) {
      p <- sys.parent(1L)
      stop(simpleError(m, call = if (p > 0L) sys.call(p)))
    }
  }
  invisible(NULL)
}

warn_if_not <- function(..., m) {
  n <- ...length()
  for (i in seq_len(n)) {
    e <- ...elt(i)
    if (!(is.logical(e) && !anyNA(e) && all(e))) {
      p <- sys.parent(1L)
      warning(simpleWarning(m, call = if (p > 0L) sys.call(p)))
    }
  }
  invisible(NULL)
}

merge_frames <- function(x) {
  if (all(vapply(x, inherits, logical(1L), "list"))) {
    d <- merge_frames(lapply(x, merge_frames))
  } else if (all(vapply(x, is.data.frame, logical(1L)))) {
    s <- if (is.null(names(x))) as.character(seq_along(x)) else names(x)
    n <- vapply(x, nrow, integer(1L))
    d <- data.frame(.X = rep(s, n), do.call(rbind, x), stringsAsFactors = TRUE)
  } else {
    stop("`x` must be a list of (lists of) data frames.")
  }
  d
}
