## Load
d <- read.csv("IPHIS_REPORT_2years.csv")
d <- d[c("COUNTY_NAME", "CASE_REPORTED_DATE")]
names(d) <- c("region", "Date")

## Clean
d$region <- factor(d$region, exclude = "")
d$Date <- as.Date(d$Date, format = "%d%b%Y:00:00:00")
d <- d[complete.cases(d), , drop = FALSE]

## Aggregate
ontario <- aggregate(seq_len(nrow(d)), d, length)
names(ontario) <- c(names(d), "cases_new")
o <- do.call(order, ontario[names(d)])
ontario <- ontario[o, , drop = FALSE]
row.names(ontario) <- NULL

## Save
save(ontario, file = "ontario.RData")
