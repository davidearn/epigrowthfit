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
#'   A named list of lists of integer, Date, or character vectors
#'   of length 2. `names(subset)` must be a subset of `levels(ts)`.
#'   `subset[[s]]` defines fitting windows for time series `s`.
#'   See Details.
#'
#' @details
#' ## Subsetting
#' In `subset[[s]]`, integer vectors must index time points in
#' time series `s`. Character and Date vectors are converted to
#' index vectors as follows: (i) all character vectors `x` are
#' replaced with `as.Date(x)`, then (ii) each Date is replaced
#' with the index of the nearest time point in time series `s`.
#'
#' The resulting set of vectors `c(i, j)` (which must satisfy
#' `i < j`) are replaced with `i:j`, and these vectors (which
#' must be disjoint) each index the set of time points included
#' in a fitting window.
#'
#' ## Interpreting `index`
#'
#' `index` has levels `seq_len(sum(lengths(subset)))`.
#' `split(date, index)` splits `date` by fitting window,
#' `is.na(index)` indexes points in `date` not belonging
#' to a window.
#'
#' `index` can be assigned to the so-named argument of [egf()].
#'
#' ## Interpreting `wave`
#'
#' `wave` has levels `seq_len(max(lengths(subset)))`.
#' It is a recoding of `index`, obtained within a
#' given time series by (i) replacing elements in
#' the `k`th fitting window of that time series
#' with `k`, (ii) replacing `NA` by carrying last
#' observations forward, then (iii) replacing any
#' remaining (leading) `NA` with 1.
#'
#' If the `k`th fitting window of each time series
#' corresponds to the `k`th epidemic wave in that
#' time series, and if epidemic wave is being included
#' as a fixed or random effect in a mixed effects model,
#' then `wave` can be passed to [egf()] as a variable
#' in `data`.
#'
#' @return
#' A data frame with `length(date)` rows,
#' listing two factors, `index` and `wave`.
#'
#' @export
make_index <- function(date, ts, subset) {
  ## Validate `date`, `ts`
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
  tts <- tabulate(ts)
  stop_if_not(
    tts >= 2L,
    m = "Level sets of `ts` must have size 2 or greater."
  )
  date_split <- split(date, ts)
  stop_if_not(
    vapply(date_split, function(x) all(diff(x) > 0), FALSE),
    m = "`date` must be increasing in each level of `ts`."
  )

  ## Validate `subset`, coercing date to index if necessary
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
  date_to_index <- function(x, d) vapply(x, function(d0) which.min(abs(d - d0)), 0L)
  for (s in names(subset)) {
    l <- subset[[s]]
    d <- date_split[[s]]
    is_character <- vapply(l, is.character, FALSE)
    l[is_character] <- lapply(l[is_character], as.Date)
    is_date <- vapply(l, inherits, FALSE, "Date")
    is_numeric <- vapply(l, is.numeric, FALSE)
    stop_if_not(
      is_date | is_numeric,
      lengths(l) == 2L,
      vapply(l, diff, 0) > 0,
      unlist(l[is_numeric]) %in% seq_along(d),
      m = paste0(
        "Elements of `subset` must list Date or integer\n",
        "index vectors (length 2, increasing, without NA)."
      )
    )
    l[is_date] <- lapply(l[is_date], date_to_index, d = d)
    stop_if_not(
      vapply(l[is_date], diff, 0) > 0,
      m = sprintf("`subset[[%s]]` contains intervals of width 0\nafter conversion from Date to index.", s)
    )
    l <- l[order(vapply(l, `[`, 0L, 1L))]
    firsts <- vapply(l, `[`, 0L, 1L)
    lasts <- mapply(`[`, l, lengths(l))
    stop_if_not(
      diff(c(rbind(firsts, lasts))) > 0,
      m = "Elements of `subset` must list disjoint intervals."
    )
    subset[[s]] <- l
  }
  subset <- subset[order(match(names(subset), levels(ts)))]

  ## Actually do stuff

  N <- length(date)
  index_levels <- seq_len(sum(lengths(subset)))
  wave_levels <- seq_len(max(lengths(subset)))
  index_split <- split(rep(factor(NA, levels = index_levels), N), ts)
  wave_split <- split(rep(factor(1, levels = wave_levels), N), ts)

  k <- 1L
  for (s in names(subset)) {
    l <- subset[[s]]
    n <- length(wave_split[[s]])

    ## Replace `c(a, b)` with `a:b`
    firsts <- vapply(l, `[`, 0L, 1L)
    lasts <- mapply(`[`, l, lengths(l))
    l <- Map(`:`, firsts, lasts)

    ## Assign levels to specified subsets of `index`, `wave`
    index_split[[s]][unlist(l)] <- rep.int(seq.int(k, length.out = length(l)), lengths(l))
    wave_split[[s]][unlist(l)] <- rep.int(seq_along(l), lengths(l))

    ## Fill in the rest of `wave` with LOCF
    i1 <- c(firsts[-1L], n + 1L)
    i2 <- lasts
    for (j in which(i1 - i2 > 1L)) {
      wave_split[[s]][seq.int(i2[j] + 1L, i1[j] - 1L)] <- j
    }
    k <- k + length(l)
  }

  data.frame(
    index = unsplit(index_split, ts),
    wave  = unsplit(wave_split, ts)
  )
}
