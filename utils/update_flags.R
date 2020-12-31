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

fn_cpp <- "../src/epigrowthfit.cpp"
fn_r   <- "../R/egf_utils_copy.R"

lines_cpp <- readLines(fn_cpp)
lines_r   <- readLines(fn_r)

curve_names <- get_enum_names(lines_cpp, "curve")
distr_names <- get_enum_names(lines_cpp, "distr")

pattern <- "^(  curve_names <- c\\().*(\\) # GREP_FLAG)$"
replacement <- sprintf("\\1%s\\2", paste(sprintf("\"%s\"", curve_names), collapse = ", "))
lines_r <- sub(pattern, replacement, lines_r)

pattern <- "^(  distr_names <- c\\().*(\\) # GREP_FLAG)$"
replacement <- sprintf("\\1%s\\2", paste(sprintf("\"%s\"", distr_names), collapse = ", "))
lines_r <- sub(pattern, replacement, lines_r)

writeLines(lines_r, fn_r)
