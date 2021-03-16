#' Safely convert decimal dates to dates
#'
#' Converts decimal dates to dates in a way that avoids truncation
#' by [as.Date()] (see Details).
#'
#' @param x A numeric vector listing decimal dates.
#'
#' @return
#' A Date vector with length equal to `length(dd)`.
#'
#' @details
#' [lubridate::date_decimal()] converts numeric vectors to POSIXct
#' objects. The naive way to convert POSIXct objects (storing both
#' date and time) to Date objects (storing just date) is to apply
#' [as.Date()] directly. However, [as.Date.POSIXct()] performs
#' this conversion by "ignoring the time after midnight in the
#' representation of the time in the specified time zone". Hence
#' the error will be undesirably large if the time is past noon.
#' The solution is to round date-times beforehand by calling
#' [round.POSIXt()] with argument `units = "days"`, before calling
#' [as.Date.POSIXct()].
#'
#' @examples
#' datetime_numeric <- 2001.002
#' datetime_POSIXct <- lubridate::date_decimal(date_numeric) # time past noon
#' datetime_POSIXct_rounded <- round(datetime_POSIXct, units = "days")
#' date_bad <- as.Date(datetime_POSIXct)
#' date_good <- as.Date(datetime_POSIXct_rounded)
#' identical(date_good, safe_Date(datetime_numeric))
#'
#' @keywords internal
safe_Date <- function(x) {
  as.Date(round(lubridate::date_decimal(x), units = "days"))
}

#' Aggregate count data
#'
#' @description
#' Aggregates count data to obtain time series of a desired spacing.
#' Uniformly distributes counts onto the dates they represent to
#' achieve equal spacing. This is useful if, for example, the count
#' data are not equally spaced or the aggregation interval is not an
#' integer multiple of the observation interval.
#'
#' @param date
#'   A Date vector listing (increasing) time points.
#' @param counts
#'   A data frame with numeric variables listing counts on the
#'   given dates, so that `nrow(counts) = length(date)`. Each
#'   variable defines a time series to be aggregated. Must be
#'   right-aligned, i.e., `counts[i, ]` must give counts over
#'   the days from `date[i-1]+1` to `date[i]`. The first row
#'   is ignored, since `date` provides no information about how
#'   it should be interpreted.
#' @param aggregation
#'   An integer indicating the number of days over which to aggregate
#'   (e.g., 7 for weekly aggregates).
#' @param alignment
#'   An element of `c("r", "c", "l")` indicating how aggregates
#'   should be aligned in the returned data frame. See Value.
#'
#' @return
#' A data frame whose first column is a Date vector listing
#' time points *after* aggregation and whose other columns
#' are aggregations of the columns in argument `counts`. The
#' dates and the correct interpretation of aggregates depend
#' on the value of argument `alignment` (`"r"`, `"c"`, or `"l"`).
#' Let `d` be the date column and `x` any other column. Then:
#'
#' \describe{
#'   \item{`"r"`}{
#'     `x[i]` is an aggregate over the days from `d[i-1]+1` to `d[i]`.
#'   }
#'   \item{`"c"`}{
#'     If `aggregation` is odd and equal to `2*m+1`, then `x[i]`
#'     is an aggregate over the days from `d[i]-m` to `d[i]+m`.
#'     If `aggregation` is even and equal to `2*n`, then `x[i]`
#'     is an aggregate over the days from `d[i]-n` to `d[i]+n-1`.
#'   }
#'   \item{`"l"`}{
#'     `x[i]` is an aggregate over the days from `d[i]` to `d[i+1]-1`.
#'   }
#' }
#'
#' @details
#' Aggregates will not be integer valued unless the aggregation
#' interval is an integer multiple of the observation interval.
#'
#' @examples
#' date <- as.Date("2001-01-01") + cumsum(sample(1:14, 100, replace = TRUE))
#' counts <- data.frame(
#'   x = rbinom(100, 12, prob = 0.2),
#'   y = rbinom(100, 12, prob = 0.5)
#' )
#' ac_out <- lapply(c("r", "c", "l"), function(x) {
#'   aggregate_counts(date, counts, aggregation = 7, alignment = x)
#' })
#' names(ac_out) <- c("r", "c", "l")
#' lapply(ac_out, head)
#'
#' @keywords internal
#' @importFrom utils head
aggregate_counts <- function(date,
                             counts,
                             aggregation = 1,
                             alignment = "r") {
  if (!is.numeric(aggregation) || length(aggregation) != 1 ||
      !isTRUE(aggregation >= 1) || aggregation %% 1 != 0) {
    stop("`aggregation` must be a positive integer.")
  }
  if (!is.character(alignment) || length(alignment) != 1 ||
      !alignment %in% c("r", "c", "l")) {
    stop("`alignment` must be \"r\", \"c\", or \"l\".")
  }
  if (!inherits(date, "Date")) {
    stop("`date` must have class \"Date\".")
  } else if (length(date) < aggregation + 1) {
    stop("`length(date)` must be greater than or equal to `aggregation + 1`.")
  } else if (!isTRUE(all(diff(date) > 0))) {
    stop("`date` must be increasing.")
  }
  if (!is.data.frame(counts) || !isTRUE(all(sapply(counts, is.numeric)))) {
    stop("`counts` must be a data frame with numeric variables.")
  } else if (isTRUE(any(counts < 0))) {
    stop("`counts` must not have negative elements.")
  } else if (nrow(counts) != length(date)) {
    stop("`nrow(counts)` must be equal to `length(date)`.")
  }

  ## Uniformly distribute counts onto aggregated days
  ## to get a daily time series
  diff_date <- as.numeric(diff(date))
  date_daily <- seq(date[2], date[length(date)], by = 1)
  counts_daily <- as.data.frame(sapply(counts[-1, , drop = FALSE], function(x) {
    rep(x / diff_date, times = diff_date)
  }))

  ## Discard unused counts
  n <- floor(length(date_daily) / aggregation) * aggregation
  date_daily <- date_daily[1:n]
  counts_daily <- head(counts_daily, n)

  ## First date counted in each aggregate
  starts <- seq(date_daily[1], date_daily[n], by = aggregation)
  ## Last date counted in each aggregate
  ends <- seq(date_daily[aggregation], date_daily[n], by = aggregation)
  ## Date assigned to each aggregate
  date <- switch(alignment,
    l = starts,
    c = starts + ceiling((ends - starts) / 2),
    r = ends
  )
  ## Aggregates
  counts <- sapply(counts_daily, function(x) {
    colSums(matrix(x, nrow = aggregation))
  })
  data.frame(date, counts)
}

#' Factor case counts by outbreak and severity
#'
#' @description
#' Assigns time points in a data frame an outbreak label,
#' then assigns outbreaks a severity label.
#'
#' @param data
#'   A data frame with a numeric variable `time` listing time points.
#' @param outbreak_definitions
#'   A data frame with variables:
#'
#'   \describe{
#'     \item{`outbreak`}{
#'       \[character or numeric\] Outbreak labels, without duplicates.
#'     }
#'     \item{`severity`}{
#'       \[factor\] Severity labels, classifying each outbreak by its
#'       severity (e.g., `"minor"` and `"major"`).
#'     }
#'     \item{`start`,`end`}{
#'       \[numeric\] Increasing vectors defining the start and
#'       end times of each outbreak. To ensure that the intervals
#'       are mutually disjoint, `c(rbind(start, end))` must also
#'       be increasing.
#'     }
#'   }
#'
#' @return
#' Argument `data` with certain rows dropped and two factor
#' variables `outbreak` and `severity` added. See Details.
#'
#' @details
#' `data$time` is factored according to membership in the intervals
#' defined by `outbreak_definitions[, c("start", "end")]`. Levels
#' for this factor are taken from `outbreak_definitions$outbreak`.
#' The resulting factor is added to `data` as `data$outbreak`, then
#' times not belonging to any interval are discarded.
#'
#' `data$outbreak` is then factored according to severity
#' using `outbreak_definitions[, c("outbreak", "severity")]`
#' as a dictionary. The resulting factor is added to `data`
#' as `data$severity`, then `data` is returned.
#'
#' @keywords internal
factor_by_outbreak <- function(data, outbreak_definitions) {
  ## Factor times by outbreak
  data$outbreak <- cut(data$time,
    breaks = c(t(outbreak_definitions[, c("start", "end")])),
    labels = c(rbind(NA, outbreak_definitions$outbreak))[-1],
    right  = FALSE,
    ordered_result = TRUE
  )
  data$outbreak <- droplevels(data$outbreak, exclude = NA)
  data <- data[!is.na(data$outbreak), ]
  row.names(data) <- NULL
  if (nrow(data) == 0) {
    stop("Union of intervals in `outbreak_definitions`, ",
         "contains 0 elements of `data$time`.")
  }

  ## Factor outbreaks by severity
  data$severity <- data$outbreak
  levels(data$severity) <- outbreak_definitions$severity[match(levels(data$outbreak), outbreak_definitions$outbreak)]
  if (all(levels(data$severity) %in% c("minor", "major"))) {
    ## Order the levels if they happen to be "minor" and "major"
    data$severity <- factor(data$severity,
      levels = c("minor", "major"),
      ordered = TRUE
    )
  }

  data
}
