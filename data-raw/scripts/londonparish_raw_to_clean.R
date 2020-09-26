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
list_of_aggregates <- by(df,
  INDICES = df$outbreak,
  FUN = function(x) {
    data <- aggregate_counts(
      date        = x$date,
      counts      = x["burials"],
      aggregation = aggregation,
      alignment   = alignment
    )
    data$time <- lubridate::decimal_date(data$date)
    data$outbreak <- x$outbreak[1]
    data$severity <- x$severity[1]
    data
  },
  simplify = FALSE
)
## Merge the aggregates (one for each level of `outbreak`)
df <- do.call(rbind, list_of_aggregates)
## Clean up
df <- df[, c("date", "time", "burials", "outbreak", "severity")]
rownames(df) <- NULL


londonparish <- structure(df, aggregation = aggregation, alignment = alignment)
save(londonparish, file = "../../data/londonparish.RData")
