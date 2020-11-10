source("utils.R")
aggregation <- 7
alignment <- "r"
load("../londonparish_raw.RData")
load("../outbreak_definitions.RData")
df <- londonparish_raw


## Factor times by outbreak and outbreaks by severity
df <- factor_by_outbreak(df, outbreak_definitions)

## In each level of `outbreak`, aggregate `wills` over
## `aggregation` days and report aggregates on the
## first, central, or last day depending on `alignment`
df_split <- split(df, df$outbreak)
df_split <- lapply(df_split, function(x) {
  y <- aggregate_counts(
    date        = x$date,
    counts      = x["burials"],
    aggregation = aggregation,
    alignment   = alignment
  )
  y$time <- lubridate::decimal_date(y$date)
  y$outbreak <- x$outbreak[1]
  y$severity <- x$severity[1]
  y
})

## Merge and clean
df <- do.call(rbind, df_split)
df <- df[, c("date", "time", "burials", "outbreak", "severity")]
row.names(df) <- NULL


londonparish <- structure(df, aggregation = aggregation, alignment = alignment)
save(londonparish, file = "../../data/londonparish.RData")
