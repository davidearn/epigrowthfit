\name{egf_parallel}
\alias{egf_parallel}
\title{Define a parallelization method}
\usage{
egf_parallel(
  method = c("serial", "multicore", "snow"),
  outfile = "",
  cores = getOption("egf.cores", 1L),
  args = list(),
  cl = NULL
)
}
\arguments{
\item{method}{A character string indicating a method of parallelization.
\code{"\link[=lapply]{serial}"} indicates no parallelization.
\code{"\link[parallel:mclapply]{multicore}"} indicates \R level forking.
It is intended for use from a terminal rather than a GUI.
On Windows, \code{"multicore"} is equivalent to \code{"serial"}.
\code{"\link[parallel:parLapply]{snow}"} indicates socket clusters.
\code{"snow"} is supported on both Unix-alikes and Windows.}

\item{outfile}{A character string indicating a file path where console output
should be diverted. An empty string indicates no diversion.
If \code{method = "snow"}, then diversion may be necessary to view output.}

\item{cores}{A positive integer indicating a number of threads/processes
to fork/spawn when \code{parallel != "serial"}. The maximum
is typically \code{\link[parallel]{detectCores}(TRUE, FALSE)}.}

\item{args}{A list of optional arguments to
\code{\link[parallel]{mclapply}} (\code{method = "multicore"}) or
\code{\link[parallel]{makePSOCKcluster}} (\code{method = "snow"}).}

\item{cl}{An existing \link[parallel:makePSOCKcluster]{socket cluster}
(\code{method = "snow"}).
The default (\code{NULL}) is to create a new clusters as necessary
and terminate them upon job completion.
(If non-\code{NULL}, then \code{outfile}, \code{cores}, and \code{args}
are ignored.)}
}
\value{
A list inheriting from class \code{"egf_parallel"}
containing the arguments (after possible matching and coercion).
}
\description{
Defines instructions for parallelization by linking a method with options.
}
\details{
For general information about parallelism in \R,
see \code{\link{vignette}("parallel", "parallel")}.
}
\examples{
parallel <- egf_parallel()
str(parallel)

}
\seealso{
Other behaviour-defining functions: 
\code{\link{egf_control}()},
\code{\link{egf_model}()},
\code{\link{egf_optimizer}()}
}
\concept{behaviour-defining functions}
