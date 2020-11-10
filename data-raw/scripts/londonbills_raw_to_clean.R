source("utils.R")
aggregation <- 7
alignment <- "r"
load("../londonbills_raw.RData")
load("../outbreak_definitions.RData")
df <- londonbills_raw


## Factor times by outbreak and outbreaks by severity
df <- factor_by_outbreak(df, outbreak_definitions)

## In each level of `outbreak`, aggregate `all_causes_deaths` and
## `plague_deaths` over `aggregation` days and report aggregates
## on the first, central, or last day depending on `alignment`
df_split <- split(df, df$outbreak)
df_split <- lapply(df_split, function(x) {
  y <- aggregate_counts(
    date        = x$date,
    counts      = x[c("all_causes_deaths", "plague_deaths")],
    aggregation = aggregation,
    alignment   = alignment
  )
  y$time <- lubridate::decimal_date(y$date)
  y$outbreak <- x$outbreak[1]
  y$severity <- x$severity[1]
  y$population <- approx(x$time, x$population, xout = y$time)$y
  y
})

## Merge and clean
df <- do.call(rbind, df_split)
df <- df[, c("date", "time", "all_causes_deaths", "plague_deaths",
             "outbreak", "severity", "population")]
row.names(df) <- NULL


londonbills <- structure(df, aggregation = aggregation, alignment = alignment)
save(londonbills, file = "../../data/londonbills.RData")
