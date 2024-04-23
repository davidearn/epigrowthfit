\name{egf_cache}
\alias{egf_cache}
\title{Cache Example Objects}
\description{
A utility for caching objects created when code in the \sQuote{Examples}
section of an \pkg{epigrowthfit} help page is sourced, typically by a
call to \code{\link{example}}.  It is not intended for use outside of
this context.  It is exported so that it can be accessed without
\code{\link{:::}} in \sQuote{Examples} and to provide a way for users
to straightforwardly clear the cache from an \R{} process.
}
\usage{
egf_cache(file, object, topic = NULL, clear = FALSE, clearAll = FALSE, \dots)
}
\arguments{
\item{file}{
  a character string containing the base name of a file.
  The absolute path to the cache file is constructed as:
\preformatted{file.path(tools::R_user_dir("epigrowthfit", "cache"),
          utils::packageVersion("epigrowthfit"),
          file)}}
\item{object}{
  an \R{} object.}
\item{topic}{
  a character string specifying the topic associated with \code{file}.
  If \code{NULL}, then an attempt is made to recover \code{topic} from
  \code{file}, and an error is thrown if the attempt fails.}
\item{clear}{
  a logical flag.  If \code{TRUE}, then \code{file} is removed by
  \code{\link{unlink}} if it exists, and nothing else is done.}
\item{clearAll}{
  a logical flag.  If \code{TRUE}, then the cache directory is
  removed by \code{\link{unlink}} if it exists, and nothing else
  is done.  \code{clearAll} is resolved before \code{clear}.}
\item{\dots}{
  optional arguments passed to \code{\link{saveRDS}}.}
}
\value{
If \code{clear = FALSE} and \code{clearAll = FALSE},
then \code{object} or an \R{} object restored from \code{file}
by \code{\link{readRDS}}.

If \code{clear = TRUE} or \code{clearAll = TRUE},
then the value of the call to \code{\link{unlink}}.
}
\section{Behaviour when \code{clear = FALSE} and \code{clearAll = FALSE}}{
If \code{file} is found by \code{\link{file.exists}}, then the
cached object is restored by \code{\link{readRDS}} and returned.
In this case, \code{object} is never evaluated.

If \code{file} is not found but \code{object} is supplied,
then \code{object} is written to \code{file} by \code{\link{saveRDS}}
and returned.

If \code{file} is not found and \code{object} is missing,
then \code{\link{example}(topic, character.only = TRUE)}
is run in order to create the file.  If this works, then
the cached object is restored using \code{\link{readRDS}}.
Otherwise, an error is thrown.
}
\section{File Naming Scheme}{
In calls to \code{egf_cache} found in the \sQuote{Examples}
for topic \code{topic}, \code{file} is set equal to the string
obtained by

\itemize{
\item replacing non-word characters in \code{topic} with hyphens,
then
\item appending a hyphen, an integer index, and a \file{.rds}
extension to the end of the resulting string.
}

This scheme allows \code{topic} to be recovered from \code{file} via
pattern matching against topics found by \code{\link{help.search}}.
}
\examples{
verdir <- file.path(tools::R_user_dir("epigrowthfit", "cache"),
                    utils::packageVersion("epigrowthfit"))
lf1 <- list.files(verdir)
file <- "test.rds"
a <- 1
x <- egf_cache(file, a)
y <- egf_cache(file, a <- 2)
lf2 <- list.files(verdir)
z <- egf_cache(file, clear = TRUE)
lf3 <- list.files(verdir)

stopifnot(identical(x, 1),
          identical(y, 1),
          identical(z, 0L),
          identical(a, 1),
          identical(setdiff(lf2, lf1), file),
          identical(lf3, lf1))
}
