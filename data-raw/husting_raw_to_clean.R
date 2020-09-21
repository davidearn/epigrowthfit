source("factor_by_outbreak.R")
source("aggregate_counts.R")
aggregation <- 7
alignment <- "r"
load("husting_raw.RData")
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


husting_daily <- df
save(husting_daily, file = "husting_daily.RData")


## Factor times by outbreak and outbreaks by severity
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
husting <- df
save(husting, file = "../data/husting.RData")
