pkg <- "epigrowthfit"
obj <- "covid19.canada"

if (!file.exists("DESCRIPTION") ||
	read.dcf("DESCRIPTION", fields = "Package")[1L] != pkg)
	stop(gettextf("script must be run in package [%s] root directory", pkg),
	     domain = NA)

url <- "https://wzmli.github.io/COVID19-Canada/COVID19_Canada.csv"
csv <- tempfile()
download.file(url, csv)
dat <- raw <- read.csv(csv)
str(raw)

stopifnot(exprs = {
	is.data.frame(dat)
	(m <- match(c("Province", "Date", "confirmed_positive"), names(dat))) > 0L
	vapply(dat[m], typeof, "") == c("character", "character", "integer")
})
dat <- dat[m]
names(dat) <- c("province", "date", "confirmed")

p <- c("AB", "BC", "MB", "NB" , "NL", "NS",
       "NT", "NU", "ON", "PEI", "QC", "SK", "YT")
dat[["province"]] <- factor(dat[["province"]],
                            levels = p,
                            labels = replace(p, p == "PEI", "PE"))
dat[["date"]] <- as.Date(dat[["date"]]) + 1 # assume right alignment

stopifnot(!anyNA(dat[c("province", "date")]))
table(dat[["province"]])

## Discard confirmations after 2022-05-01 00:00:00
s <- dat[["date"]] <= as.Date("2022-05-01")
dat <- dat[s, , drop = FALSE]

## Order by province, then date, then confirmed
o <- do.call(order, unname(dat))
dat <- dat[o, , drop = FALSE]
is.unsorted(o)

tapply(dat[["date"]], dat[["province"]],
       function(x) sum(duplicated(x)))
tapply(dat[["confirmed"]], dat[["province"]],
       function(x) sum(is.na(x)))

f <- function(d) {
	## Discard leading NA and duplicated times (keeping argmax)
	n <- nrow(d)
	k <- seq.int(which.min(is.na(d[["confirmed"]])), n)
	k <- k[!duplicated(d[["date"]][k], fromLast = TRUE)]
	d[k, , drop = FALSE]
}
dat <- do.call(rbind, by(dat, dat[["province"]], f))
row.names(dat) <- NULL

assign(obj, dat)
save(list = obj, file = file.path("data", paste0(obj, ".rda")))
str(dat)
