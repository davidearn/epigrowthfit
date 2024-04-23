\name{egf_parallel}
\alias{egf_parallel}
\title{Define a Parallelization Method}
\description{
Defines instructions for parallelization by linking a method with options.
}
\usage{
egf_parallel(method = c("serial", "multicore", "snow"),
             outfile = "", cores = getOption("egf.cores", 1L),
             args = list(), cl = NULL)
}
\arguments{
\item{method}{
  a character string indicating a method of parallelization.
  \code{"serial"} indicates no parallelization.
  \code{"multicore"} indicates \R{} level forking.
  It is intended for use from a terminal rather than from a GUI.
  \code{"snow"} indicates socket clusters.
  On Windows, \code{"multicore"} is equivalent to \code{"serial"}.
  \code{"snow"} is supported on both Unix-alikes and Windows.}
\item{outfile}{
  a character string indicating a file path where console output
  should be diverted.  An empty string indicates no diversion.
  If \code{method = "snow"}, then diversion may be necessary to
  view output.}
\item{cores}{
  a positive integer indicating a number of threads/processes to
  fork/spawn when \code{parallel != "serial"}.
  \code{\link[parallel]{detectCores}} can be called to detect the
  theoretical maximum.}
\item{args}{
  a list of optional arguments to
  \code{\link[parallel]{mclapply}} (\code{method = "multicore"}) or
  \code{\link[parallel]{makePSOCKcluster}} (\code{method = "snow"}).}
\item{cl}{
  an existing \link[parallel:makePSOCKcluster]{socket cluster}
  (\code{method = "snow"}).  The default is to create a new cluster
  stop it upon job completion.}
}
\value{
A list inheriting from class \code{"egf_parallel"}
containing the arguments (after possible matching and coercion).
}
\seealso{
\code{\link{vignette}("parallel", "parallel")}.
}
\examples{
parallel <- egf_parallel()
str(parallel)
}
