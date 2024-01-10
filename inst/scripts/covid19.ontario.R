pkg <- "epigrowthfit"
obj <- "covid19.ontario"

if (!file.exists("DESCRIPTION") ||
	read.dcf("DESCRIPTION", fields = "Package")[1L] != pkg)
	stop(gettextf("script must be run in package [%s] root directory", pkg),
	     domain = NA)

## Otherwise 403
options(HTTPUserAgent = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:121.0) Gecko/20100101 Firefox/121.0")

## From https://data.ontario.ca/dataset/status-of-covid-19-cases-in-ontario
url <- "https://data.ontario.ca/datastore/dump/ed270bb8-340b-41f9-a7c6-e8ef587e6d11?bom=True"
csv <- tempfile()
download.file(url, csv)

dat <- raw <- read.csv(csv, check.names = FALSE)
str(raw)

nms <- c("Reported Date", "Total Cases",
         "Total patients approved for testing as of Reporting Date",
         "Total tests completed in the last day")
stopifnot(exprs = {
	is.data.frame(dat)
	match(nms, names(dat), 0L) > 0L
	vapply(dat[nms], typeof, "") == c("character", rep.int("integer", 3L))
})

dat <- dat[nms]
names(dat) <- c("date", "cases", "00tests", "tests")
dat[["date"]] <- as.Date(dat[["date"]], format = "%Y-%m-%dT00:00:00")

stopifnot(!anyNA(dat[["date"]]))

## Discard counts after 2022-05-01 00:00:00
s <- dat[["date"]] <= as.Date("2022-05-01")
dat <- dat[s, , drop = FALSE]

stopifnot(exprs = {
	identical(range(dat[["date"]]), as.Date(c("2020-02-06", "2022-05-01")))
	identical(range(diff(dat[["date"]])), .difftime(c(1, 1), "days"))
})

## Account for changes in reporting of tests on two dates;
## see ./covid19.ontario.definitions.tsv
i0 <- match("2020-03-30", dat[["date"]], 0L)
i1 <- match("2020-04-15", dat[["date"]], 0L)
i2 <- nrow(dat) + 1L
dat[["tests"]][1L:(i0 - 1L)] <- NA
dat[["tests"]][i0:(i1 - 1L)] <- dat[["00tests"]][i0:(i1 - 1L)]
dat[["tests"]][i1:(i2 - 1L)] <- cumsum(dat[["tests"]][i1:(i2 - 1L)]) + dat[["tests"]][i1 - 1L]
dat[["00tests"]] <- NULL

stopifnot(exprs = {
	identical(which(is.na(dat[["cases"]])), 1L:2L)
	identical(which(is.na(dat[["tests"]])), 1L:(i0 - 1L))
})

## Discard leading NA
dat <- dat[-(1L:min(2L, i0 - 1L)), , drop = FALSE]

## Force automatic row names
row.names(dat) <- NULL

assign(obj, dat)
save(list = obj, file = file.path("data", paste0(obj, ".rda")))
str(dat)
