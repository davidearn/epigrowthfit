filename_h <- "../src/enum.h"
filename_r   <- "../R/enum.R"
lines_h <- readLines(filename_h)
lines_r   <- readLines(filename_r)

## Function to extract enumerator names from lines
## of file declaring enumerated type
get_enum_names <- function(lines, type) {
  start <- sprintf("^enum %s", type)
  i_start <- grep(start, lines)[1L]
  i_end <- i_start
  while (!grepl("};$", lines[i_end])) {
    i_end <- i_end + 1L
  }
  text <- paste(lines[i_start:i_end], collapse = "")
  text <- sub(start, "", text)
  text <- gsub("[{[:space:]};]", "", text) # leaves comma-separated list
  strsplit(text, ",")[[1L]]
}

## Character vectors listing enumerator names
curve_names <- get_enum_names(lines_h, "curve")
distr_names <- get_enum_names(lines_h, "distr")

## String giving quote-delimited, comma-separated list
curve_str <- paste(sprintf("\"%s\"", curve_names), collapse = ", ")
distr_str <- paste(sprintf("\"%s\"", distr_names), collapse = ", ")

## Replace argument of c() in matched line with string
lines_r <- sub(
  pattern = "^(  curve_names <- c\\().*(\\) # GREP_FLAG)$",
  replacement = sprintf("\\1%s\\2", curve_str),
  x = lines_r
)
lines_r <- sub(
  pattern = "^(  distr_names <- c\\().*(\\) # GREP_FLAG)$",
  replacement = sprintf("\\1%s\\2", distr_str),
  x = lines_r
)
writeLines(lines_r, filename_r)
