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
#' @return
#' A data frame with `length(date)` rows,
#' listing two factors, `index` and `wave`.
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
  wave_split <- split(rep(factor(1, levels = seq_len(max(lengths(subset)))), length(date)), ts)

  k <- 1L
  for (s in names(subset)) {
    l <- subset[[s]]
    stop_if_not(
      vapply(l, is.numeric, FALSE),
      lengths(l) > 0L,
      unlist(l) %in% seq_len(tts[[s]]),
      !duplicated(unlist(l)),
      m = sprintf("`subset[[%s]]` must list disjoint subsets\nof `seq_len(table(ts)[[%s]])`.", s, s)
    )
    l <- lapply(l, sort.int)
    stop_if_not(
      vapply(l, function(x) all(diff(x) == 1L), FALSE),
      m = sprintf("Elements of `subset[[%s]]` must be contiguous integer vectors.", s)
    )

    l <- l[order(vapply(l, `[[`, 0L, 1L))]
    index_split[[s]][unlist(l)] <- rep.int(seq.int(k, length.out = length(l)), lengths(l))
    wave_split[[s]][unlist(l)] <- rep.int(seq_along(l), lengths(l))

    i1 <- c(vapply(l, `[`, 0L, 1L), length(wave_split[[s]]) + 1L)
    i2 <- c(0L, vapply(l, function(i) i[length(i)], 0L))
    w <- c(1L, seq_along(l))
    for (j in seq_along(w)) {
      if (i1[j] - i2[j] > 1L) {
        wave_split[[s]][seq.int(i2[j] + 1L, i1[j] - 1L)] <- w[j]
      }
    }
    k <- k + length(l)
  }

  data.frame(
    index = unsplit(index_split, ts),
    wave  = unsplit(wave_split, ts)
  )
}
