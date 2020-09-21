source("factor_by_outbreak.R")
source("aggregate_counts.R")
aggregation <- 7
alignment <- "r"
load("london_parish_raw.RData")
df <- london_parish_raw


## Factor times by outbreak and outbreaks by severity
df <- factor_by_outbreak(df,
  outbreak_definitions = local({load("outbreak_definitions.RData"); outbreak_definitions})
)
## In each level of `outbreak`, aggregate `burials` over
## `aggregation` days and report aggregates on the
## first, central, or last day depending on `alignment`
df_by_outbreak <- by(df,
  INDICES = df$outbreak,
  FUN = function(x) {
    y <- aggregate_counts(
      date        = x$date,
      count       = x["burials"],
      aggregation = aggregation,
      alignment   = alignment
    )
    y$time <- lubridate::decimal_date(y$date)
    y$outbreak <- x$outbreak[1]
    y$severity <- x$severity[1]
    y[, names(x)]
  },
  simplify = FALSE
)
df <- do.call(rbind, df_by_outbreak)
rownames(df) <- NULL


londonparish <- structure(df, aggregation = aggregation, alignment = alignment)
save(londonparish, file = "../data/londonparish.RData")
