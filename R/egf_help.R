#' Index fitting windows
#'
#' Constructs a factor splitting time series in long format
#' by fitting window, to be passed to [egf()].
#'
#' @param date
#'   A Date vector listing time points from one or more time series.
#' @param ts
#'   A factor of length `length(date)` such that `split(date, ts)`
#'   splits `date` by time series.
#' @param subset
#'   A named list of lists of integer vectors. `names(subset)`
#'   must be a subset of `levels(ts)`. Integer vectors listed
#'   in `subset[[name]]` must index disjoint segments of the
#'   corresponding time series. That is, they must be disjoint
#'   subsets of `seq_len(table(ts)[[name]])`.
#'
#' @details
#' `index` has levels `seq_len(sum(lengths(subset)))`.
#' `split(date, index)` splits `date` by fitting window,
#' `is.na(index)` indexes points in `date` not belonging
#' to a window.
#'
#' `index` can be assigned to the so-named argument of [egf()].
#'
#' `wave` has levels `c(0, seq_len(max(lengths(subset)))`.
#' It is a recoding of `index`, replacing `NA` with `0`
#' and replacing elements in the `i`th fitting window of
#' each time series (chronologically) with `i`.
#'
#' If the `i`th fitting window of each time series corresponds
#' to the `i`th epidemic wave in that time series, and if
#' epidemic wave is being included as a fixed or random effect
#' in a mixed effects model, then `wave` can be passed to
#' [egf()] as a variable in `data`.
#'
#' Note that `wave` codes points not belonging to a fitting window
#' as `0`, regardless of epidemic wave. Since they do not affect
#' the fitted model, their coding is arbitrary, though `NA` is not
#' used as `NA` are not tolerated by [egf()] except in the incidence
#' variable.
#'
#' @return
#' A data frame listing 2 factors of length `length(date)`,
#' `index` and `wave`.
#'
#' @export
make_index <- function(date, ts, subset) {
  stop_if_not(
    inherits(date, "Date"),
    m = "`date` must be a Date vector."
  )
  stop_if_not(
    is.factor(ts),
    m = "`ts` must be a factor."
  )
  stop_if_not(
    length(date) > 0L,
    length(ts) == length(date),
    m = "`date` and `ts` must have equal, nonzero lengths."
  )
  ts <- droplevels(ts, exclude = NA)
  stop_if_not(
    !anyNA(date),
    !anyNA(ts),
    m = "`date` and `ts` must not have missing values."
  )
  stop_if_not(
    vapply(split(date, ts), function(x) all(diff(x) > 0), FALSE),
    m = "`date` must be increasing in each level of `ts`."
  )
  stop_if_not(
    is.list(subset),
    length(subset) > 0L,
    !is.null(names(subset)),
    names(subset) %in% levels(ts),
    !duplicated(names(subset)),
    vapply(subset, is.list, FALSE),
    m = paste0(
      "`subset` must be a named list of lists,\n",
      "with `names(subset)` a subset of `levels(ts)`."
    )
  )

  tts <- table(ts)
  subset <- subset[order(match(names(subset), levels(ts)))]
  index_split <- split(rep(factor(NA, levels = seq_len(sum(lengths(subset)))), length(date)), ts)
  wave_split <- split(rep(factor(0L, levels = c(0L, seq_len(max(lengths(subset))))), length(date)), ts)

  i <- 1L
  for (s in names(subset)) {
    l <- subset[[s]]
    stop_if_not(
      vapply(l, is.numeric, FALSE),
      lengths(l) > 0L,
      unlist(l) %in% seq_len(tts[[s]]),
      !duplicated(unlist(l)),
      m = sprintf("`subset[[%s]]` must list disjoint subsets\nof `seq_len(table(ts)[[%s]])`.", s, s)
    )
    stop_if_not(
      vapply(l, function(x) all(diff(sort.int(x)) == 1L), FALSE),
      m = sprintf("Elements of `subset[[%s]]` must be contiguous integer vectors.", s)
    )
    l <- l[order(vapply(l, `[[`, 1L))]
    index_split[[s]][unlist(l)] <- rep.int(seq.int(i, length.out = length(l)), lengths(l))
    wave_split[[s]][unlist(l)] <- rep.int(seq_along(l), lengths(l))
    i <- i + length(l)
  }

  data.frame(
    index = unsplit(index_split, ts),
    wave  = unsplit(wave_split, ts)
  )
}
