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
list_of_aggregates <- by(df,
  INDICES = df$outbreak,
  FUN = function(x) {
    data <- aggregate_counts(
      date        = x$date,
      counts      = x[c("all_causes_deaths", "plague_deaths")],
      aggregation = aggregation,
      alignment   = alignment
    )
    data$time <- lubridate::decimal_date(data$date)
    data$outbreak <- x$outbreak[1]
    data$severity <- x$severity[1]
    data$population <- approx(x$time, x$population, xout = data$time)$y
    data
  },
  simplify = FALSE
)
## Merge the aggregates (one for each level of `outbreak`)
df <- do.call(rbind, list_of_aggregates)
## Clean up
df <- df[, c("date", "time", "all_causes_deaths", "plague_deaths",
             "outbreak", "severity", "population")]
rownames(df) <- NULL


londonbills <- structure(df, aggregation = aggregation, alignment = alignment)
save(londonbills, file = "../../data/londonbills.RData")
