source("utils.R")
aggregation <- 7
alignment <- "r"
load("../husting_raw.RData")
load("../outbreak_definitions.RData")
df <- husting_raw


## Keep write dates and ignore probate dates
df <- df[, grepl("written", names(df))]

## Omit all 46 wills with a missing date element (year, month, or day)
sum(apply(df, 1, anyNA)) # 46
df <- na.omit(df)

## Count number of wills on each date
counts <- table(df$date_written)
df <- data.frame(
  date  = as.Date(dimnames(counts)[[1]]),
  wills = as.numeric(counts)
)

## Fill with zeros on missing dates
dates_all <- with(df, seq(min(date), max(date), by = 1))
dates_missing <- dates_all[!dates_all %in% df$date]
df <- rbind(df, data.frame(date = dates_missing, wills = 0))
df <- df[order(df$date), ]

## Factor times by outbreak and outbreaks by severity
df$time <- lubridate::decimal_date(df$date)
df <- factor_by_outbreak(df, outbreak_definitions)

## In each level of `outbreak`, aggregate `wills` over
## `aggregation` days and report aggregates on the
## first, central, or last day depending on `alignment`
df_split <- split(df, df$outbreak)
df_split <- lapply(df_split, function(x) {
  y <- aggregate_counts(
    date        = x$date,
    counts      = x["wills"],
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
df <- df[, c("date", "time", "wills", "outbreak", "severity")]
row.names(df) <- NULL


husting <- structure(df, aggregation = aggregation, alignment = alignment)
save(husting, file = "../../data/husting.RData")
