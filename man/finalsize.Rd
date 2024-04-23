\name{finalsize}
\alias{finalsize}
\title{Compute the Expected Epidemic Final Size}
\description{
Computes the proportion of a population expected to be infected over
the course of an epidemic, as a function of the basic reproduction number.
}
\usage{
finalsize(R0, S0, I0)
}
\arguments{
\item{R0}{
  a numeric vector listing non-negative values for the basic reproduction
  number.}
\item{S0, I0}{
  numeric vectors listing values in the interval \eqn{[0,1]} for the
  proportions of the population that are susceptible and infected,
  respectively, at the start of the epidemic.  Hence \code{S0 + I0}
  must be less than or equal to 1.}
}
\value{
A numeric vector listing values in the interval \eqn{[0,1]} for the
expected epidemic final size.
}
\details{
At least one of \code{S0} and \code{I0} must be supplied.
If \code{S0} (\code{I0}) is supplied but not \code{I0} (\code{S0}),
then the latter is assigned the value of one minus the former.

\code{R0}, \code{S0}, and \code{I0} are recycled to a common length
(the maximum of their lengths).
}
\section{Computation}{
The basic reproduction number \code{R0} defines the expected
epidemic final size \code{Z} through an implicit equation,
\preformatted{Z = S0 * (1 - exp(-R0 * (Z + I0))) ,}
which admits an explicit solution
\preformatted{Z = S0 + (1/R0) * W(-R0 * S0 * exp(-R0 * (S0 + I0))) .}
Here, \code{W} denotes the
\href{https://en.wikipedia.org/wiki/Lambert_W_function}{
  Lambert \eqn{W} function}.
\code{finalsize} computes this solution, relying on function
\code{lambertW} from package \pkg{emdbook}.
}
\seealso{
\code{\link{timescale}}, \code{\link{R0}}.
}
\references{
Ma, J. & Earn, D. J. D. (2006).
Generality of the final size formula for an epidemic of a newly invading
infectious disease.
\emph{Bulletin of Mathetmatical Biology},
\emph{68}(3), 679-702.
\doi{10.1007/s11538-005-9047-7}
}
\examples{
R0 <- 10^seq(-3, 1, length.out = 151L)
plot(R0, finalsize(R0, S0 = 1, I0 = 0), type = "l", las = 1,
     xlab = "basic reproduction number",
     ylab = "final size")
}
