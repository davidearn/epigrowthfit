load("outbreak_definitions.RData")
outbreak <- outbreak_definitions$outbreak
severity <- outbreak_definitions$severity
breaks <- .Date(unlist(outbreak_definitions[c("start", "end")], FALSE, FALSE))
labels <- rep_len(NA_character_, length(breaks) - 1L)
labels[seq.int(1L, length(labels), by = 2L)] <- as.character(outbreak)
f <- function(d) {
  d$outbreak <- droplevels(cut(d$Date, breaks = breaks, labels = labels, ordered_result = TRUE), exclude = NA)
  d$severity <- severity[match(d$outbreak, outbreak)]
  d
}

load("canterbury_raw.RData")
canterbury <- aggregate(seq_len(nrow(canterbury_raw)), canterbury_raw["Date"], length)
names(canterbury) <- c("Date", "wills")
o <- order(canterbury$Date)
canterbury <- f(canterbury[o, , drop = FALSE])
row.names(canterbury) <- NULL
save(canterbury, file = "../data/canterbury.RData", compress = "xz")

load("husting_raw.RData") # use write dates, ignore probate dates
husting_raw <- husting_raw[c("year_written", "month_written", "day_written")]
names(husting_raw) <- sub("_written$", "", names(husting_raw))
husting_raw$Date <- as.Date(do.call(ISOdate, husting_raw))
husting <- aggregate(seq_len(nrow(husting_raw)), husting_raw["Date"], length)
names(husting) <- c("Date", "wills")
o <- order(husting$Date)
husting <- f(husting[o, , drop = FALSE])
row.names(husting) <- NULL
save(husting, file = "../data/husting.RData", compress = "xz")

load("londonparishes_raw.RData")
londonparishes <- f(londonparishes_raw)
save(londonparishes, file = "../data/londonparishes.RData", compress = "xz")

load("londonbills_raw.RData")
londonbills <- f(londonbills_raw)
save(londonbills, file = "../data/londonbills.RData", compress = "xz")
