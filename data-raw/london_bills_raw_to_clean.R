source("factor_by_outbreak.R")
source("aggregate_counts.R")
aggregation <- 7
alignment <- "r"
load("london_bills_raw.RData")
df <- london_bills_raw


## Factor times by outbreak and outbreaks by severity
df <- factor_by_outbreak(df,
  outbreak_definitions = local({load("outbreak_definitions.RData"); outbreak_definitions})
)
## In each level of `outbreak`, aggregate `all_causes_deaths` and
## `plague_deaths` over `aggregation` days and report aggregates
## on the first, central, or last day depending on `alignment`
df_by_outbreak <- by(df,
  INDICES = df$outbreak,
  FUN = function(x) {
    y <- aggregate_counts(
      date        = x$date,
      count       = x[c("all_causes_deaths", "plague_deaths")],
      aggregation = aggregation,
      alignment   = alignment
    )
    y$time <- lubridate::decimal_date(y$date)
    y$outbreak <- x$outbreak[1]
    y$severity <- x$severity[1]
    y$population <- approx(x$time, x$population, xout = y$time)$y
    y[, names(x)]
  },
  simplify = FALSE
)
df <- do.call(rbind, df_by_outbreak)
rownames(df) <- NULL


londonbills <- structure(df, aggregation = aggregation, alignment = alignment)
save(londonbills, file = "../data/londonbills.RData")
