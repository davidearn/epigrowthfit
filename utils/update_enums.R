h <- c("../src/enums.h", "../tests/testthat/enums.h")
r <- "../R/enums.R"

x <- paste(unlist(lapply(h, readLines)), collapse = "\n")
x <- strsplit(gsub("[ \t\n]", "", x), "enum")[[1L]][-1L]

a <- vapply(x, gsub, "", pattern = "^.*\\{(.*)\\};$", replacement = "\\1")
names(a) <- vapply(x, gsub, "", pattern = "^(.*)\\{.*$", replacement = "\\1")

f <- function(x) sprintf("c(%s)", paste(dQuote(strsplit(x, ",")[[1L]], FALSE), collapse = ", "))
cc <- vapply(a, f, "")

y <- paste(readLines(r), collapse = "\n")
for (s in names(cc)) {
  y <- sub(
    pattern = sprintf("^(.*%s[ ]*=[ ]*)c\\(.*?\\)(.*)$", s),
    replacement = sprintf("\\1%s\\2", cc[[s]]),
    y
  )
}
writeLines(y, r)

