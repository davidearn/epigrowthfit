aggregate_counts <- function(date,
                             count,
                             first = as.numeric(date[2] - date[1]),
                             aggregation = 1,
                             alignment = "r") {
  if (!is.numeric(aggregation) || length(aggregation) != 1 ||
        !isTRUE(aggregation >= 1)) {
    stop("`aggregation` must be a numeric scalar and greater than or equal to 1.")
  }
  if (missing(date)) {
    stop("Missing argument `date`.")
  } else if (!inherits(date, "Date")) {
    stop("`date` must have class \"Date\".")
  } else if (length(date) < floor(aggregation)) {
    stop("`date` must have length greater than or equal to ",
         "`floor(aggregation)`.")
  } else if (!isTRUE(all(diff(date) > 0))) {
    stop("`date` must be increasing.")
  }
  if (missing(count)) {
    stop("Missing argument `count`.")
  } else if (!(is.data.frame(count) && isTRUE(all(sapply(count, is.numeric))))) {
    stop("`count` must be numeric or a data frame with numeric variables.")
  } else if (isTRUE(any(count < 0))) {
    stop("`count` must not have negative elements.")
  } else if (nrow(count) != length(date)) {
    stop("`nrow(count)` must be equal to `length(date)`.")
  }
  if (!is.character(alignment) || length(alignment) != 1 ||
        !alignment %in% c("l", "c", "r")) {
    stop("`alignment` must be \"l\", \"c\", or \"r\".")
  }
  if (!is.numeric(first) || length(first) != 1 || !isTRUE(first >= 1)) {
    stop("`first` must be a numeric scalar and greater than or equal to 1.")
  }

  ## Uniformly distribute counts onto aggregated days
  ## to get a daily time series
  diff_date <- c(floor(first), as.numeric(diff(date)))
  date <- seq(date[1], date[length(date)], by = 1)
  count <- as.data.frame(sapply(count, function(x) rep(x / diff_date, times = diff_date)))

  ## Need to aggregate over an integer number of days
  aggregation <- floor(aggregation)

  ## Discard unused counts
  n <- floor(length(date) / aggregation) * aggregation
  date <- date[1:n]
  count <- head(count, n)

  ## First date counted in each aggregate
  starts <- seq(date[1], date[n], by = aggregation)
  ## Last date counted in each aggregate
  ends <- seq(date[aggregation], date[n], by = aggregation)
  ## Date assigned to each aggregate
  date <- switch(alignment,
    l = starts,
    c = starts + ceiling((ends - starts) / 2),
    r = ends
  )
  ## Aggregates
  count <- sapply(count, function(x) colSums(matrix(x, nrow = aggregation)))
  data.frame(date, count)
}
