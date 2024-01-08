\name{plot.egf_profile}
\alias{plot.egf_profile}
\title{Plot likelihood profiles}
\usage{
\method{plot}{egf_profile}(x, level = attr(x, "level"), sqrt = FALSE, subset = NULL, ...)
}
\arguments{
\item{x}{An \code{"\link[=profile.egf]{egf_profile}"} object.}

\item{level}{A numeric vector with elements in (0,1). If \code{sqrt = FALSE},
then line segments are drawn to show the intersection of the profile
with lines at \code{deviance = \link{qchisq}(level, df = 1)}.}

\item{sqrt}{A logical flag.
If \code{TRUE}, then square root-transformed deviance is plotted.}

\item{subset}{An expression to be evaluated in \code{x}.
It must evaluate to a valid index vector for the rows of \code{x}.
Only indexed profiles are plotted.
The default (\code{NULL}) is to plot all profiles.}

\item{...}{Graphical parameters passed to \code{\link{plot}}.}
}
\value{
\code{NULL} (invisibly).
}
\description{
A method for plotting likelihood profiles.
}
\details{
See topic \code{\link{egf_eval}} for details on nonstandard evaluation
of \code{subset}.
}
\examples{
x <- egf_cache("profile-egf-1.rds")
plot(x, type = "o", bty = "u", las = 1, main = "")

}
