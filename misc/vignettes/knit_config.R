## ---- knit_config

## Tell knitr to insert its preamble after \input{preamble.tex},
## since it won't find \documentclass{<stuff>} in the .Rnw
knitr::knit_patterns$set(header.begin = "(^|\n)[^%]*\\s*\\\\input\\{preamble.tex\\}")

## Make some substitutions in the .tex output in order to:
## * Prevent lineno from messing up breaking of chunks over pages.
## * Prevent automatic indentation after chunks.
## * Dispense with compile errors due to knitr-xcolor interaction.
##   See https://tex.stackexchange.com/questions/148188/.
knitr::knit_hooks$set(document = function(x) {
  x <- gsub("\\begin{knitrout}", "\\nolinenumbers\\begin{knitrout}", x, fixed = TRUE)
  x <- gsub("\\end{knitrout}", "\\end{knitrout}\\linenumbers\\noindent", x, fixed = TRUE)
  gsub("\\usepackage[]{color}", "\\usepackage{xcolor}", x, fixed = TRUE)
})
## Set number of digits printed in chunk output with chunk option `digits`
knitr::knit_hooks$set(digits = function(before, options, envir) {
  if (before) {
    options(digits = options$digits)
  }
})
## Use custom palette for syntax highlighting
knitr::knit_theme$set(knitr::knit_theme$get("knit_theme.css"))
## Set default values of chunk options
knitr::opts_chunk$set(
  ## if FALSE, chunk is evaluated from scratch with compile
  cache = TRUE,
  ## if FALSE, chunk is not displayed
  echo = TRUE,
  ## if FALSE, chunk is not evaluated
  eval = TRUE,
  ## if FALSE, chunk output is not displayed and `error <- FALSE`
  include = TRUE,
  ## if FALSE, evaluation stops on errors
  error = TRUE,
  ## if FALSE, warnings printed in console, not document
  warning = TRUE,
  ## if FALSE, messages printed in console, not document
  message = TRUE,
  ## number of digits printed in chunk output
  digits = 7,
  ## column width of chunk text
  column.width = 74,
  ## figure position in document
  fig.pos = "H",
  ## figure alignment on page
  fig.align = "center",
  ## graphics device
  dev = "pdf",
  ## base point size in graphics device
  dev.args = list(pointsize = 10),
  ## chunk options to be evaluated after chunk text
  eval.after = "fig.cap",
  ## if TRUE, white space around chunks is reduced
  strip.white = TRUE
)
## Set values of global options
options(
  ## column width of chunk line
  width = 70,
  ## directional quotes
  useFancyQuotes = FALSE
)
