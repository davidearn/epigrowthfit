\name{egf_cache}
\alias{egf_cache}
\title{Cache example objects}
\usage{
egf_cache(file, object, topic = NULL, clear = FALSE, clear_all = FALSE, ...)
}
\arguments{
\item{file}{A character string containing the base name of a file.
The absolute path to the cache file is
\code{\link{file.path}(tools::\link{R_user_dir}("epigrowthfit", "cache"), file)}.}

\item{object}{An \R object.}

\item{topic}{A character string specifying the topic associated with \code{file}.
If \code{NULL}, then an attempt is made to recover \code{topic}
from \code{file}, and an error is thrown if the attempt fails.}

\item{clear}{A logical flag. If \code{TRUE}, then \code{file} is deleted by
\code{\link{file.remove}} if it exists, and nothing else is done.}

\item{clear_all}{A logical flag. If \code{TRUE}, then the package subdirectory
containing all cache files is deleted if it exists, and nothing
else is done. \code{clear_all} is resolved before \code{clear}.}

\item{...}{Arguments passed to \code{\link{saveRDS}}.}
}
\value{
If \code{clear = FALSE} and \code{clear_all = FALSE},
then \code{object} or an \R object restored from \code{file}
by \code{\link{readRDS}}.

If \code{clear = TRUE} or \code{clear_all = TRUE},
then the result of \code{\link{unlink}}.
}
\description{
A utility for caching objects created when code in the Examples section
of an \pkg{epigrowthfit} help page is sourced by \code{\link{example}}.
It is not intended for use outside of this context.
(It is only exported so that it can be accessed without \code{\link{:::}}
in Examples.)
}
\section{Behaviour when clear = FALSE and clear_all = FALSE}{

If \code{file} is found by \code{\link{file.exists}}, then the
cached object is restored by \code{\link{readRDS}} and returned.
Due to lazy evaluation of arguments by \R, \code{object} remains
a \link{promise} and is never evaluated.

If \code{file} is not found but \code{object} is supplied,
then \code{object} is written to \code{file} by \code{\link{saveRDS}}
and returned.

If \code{file} is not found and \code{object} is missing,
then \code{\link{example}(topic, character.only = TRUE)}
is run in order to create the file.
If this works, then the cached object is restored using
\code{\link{readRDS}}. Otherwise, an error is thrown.
}

\section{File naming system}{

In calls to \code{egf_cache} found in the Examples for topic \code{topic}
(except for this one) \code{file} is set equal to the string obtained by
(i) replacing non-word characters in \code{topic} with hyphens; then
(ii) appending a hyphen, an integer index, and an \code{".rds"} extension
to the end of the resulting string.
This system allows \code{topic} to be recovered from \code{file} via
pattern matching against topics found by \code{\link{help.search}}.
}

\examples{
subdir <- tools::R_user_dir("epigrowthfit", "cache")
lf1 <- list.files(subdir)
file <- "test.rds"
a <- 1
x <- egf_cache(file, a)
y <- egf_cache(file, a <- 2)
lf2 <- list.files(subdir)
z <- egf_cache(file, clear = TRUE)
lf3 <- list.files(subdir)

stopifnot(identical(x, 1),
          identical(y, 1),
          identical(z, 0L),
          identical(a, 1),
          identical(setdiff(lf2, lf1), file),
          identical(lf3, lf1))

}
