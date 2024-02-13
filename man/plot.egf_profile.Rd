\name{plot.profile.egf}
\alias{plot.profile.egf}
\title{Plot likelihood profiles}
\usage{
\method{plot}{profile.egf}(x, level = attr(x, "level"), sqrt = FALSE, subset = NULL, ...)
}
\arguments{
\item{x}{An \code{"\link{profile.egf}"} object.}

\item{level}{A numeric vector with elements in (0,1). If \code{sqrt = FALSE},
then line segments are drawn to show the intersection of the profile
with lines at \code{deviance = \link{qchisq}(level, df = 1)}.}

\item{sqrt}{A logical flag.
If \code{TRUE}, then square root-transformed deviance is plotted.}

\item{subset}{An index vector for \code{seq_len(nrow(x))}
or a language object evaluating in \code{x} to such a vector.
\code{subset} indicates which profiles are plotted;
the default (\code{NULL}) indicates all of them.}

\item{...}{Graphical parameters passed to \code{\link{plot}}.}
}
\value{
\code{NULL} (invisibly).
}
\description{
A method for plotting likelihood profiles.
}
\examples{
x <- egf_cache("profile-egf-1.rds")
plot(x, type = "o", bty = "u", las = 1, main = "")

}
