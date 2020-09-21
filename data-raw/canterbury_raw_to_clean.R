load("canterbury_raw.RData")
df <- canterbury_raw
source("factor_by_outbreak.R")
source("aggregate_counts.R")
aggregation <- 7
alignment <- "r"


## Count number of wills on each date
counts <- table(df$date)
df <- data.frame(
  date  = as.Date(dimnames(counts)[[1]]),
  wills = as.vector(counts)
)
## Fill with zeros on missing dates
dates_all <- with(df, seq(min(date), max(date), by = 1))
dates_missing <- dates_all[!dates_all %in% df$date]
df <- rbind(df, data.frame(
  date  = dates_missing,
  wills = 0
))
df <- df[order(df$date), ]
## Clean up
df$time <- lubridate::decimal_date(df$date)
df <- df[, c("date", "time", "wills")]
rownames(df) <- NULL


canterbury_daily <- df
save(canterbury_daily, file = "canterbury_daily.RData")


## Factor times by outbreak and outbreaks by severity
df$time <- lubridate::decimal_date(df$date)
df <- factor_by_outbreak(df,
  outbreak_definitions = local({load("outbreak_definitions.RData"); outbreak_definitions})
)
## In each level of `outbreak`, aggregate `wills` over
## `aggregation` days and report aggregates on the
## first, central, or last day depending on `alignment`
df_by_outbreak <- by(df,
  INDICES = df$outbreak,
  FUN = function(x) {
    y <- aggregate_counts(
      date        = x$date,
      count       = x["wills"],
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


df <- structure(df, aggregation = aggregation, alignment = alignment)
canterbury <- df
save(canterbury, file = "../data/canterbury.RData")
