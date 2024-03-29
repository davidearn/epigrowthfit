\name{finalsize}
\alias{finalsize}
\title{Compute the expected epidemic final size}
\usage{
finalsize(R0, S0, I0)
}
\arguments{
\item{R0}{A numeric vector listing non-negative values for the basic reproduction
number.}

\item{S0, I0}{Numeric vectors listing values in the interval [0,1] for the proportions
of the population that are susceptible and infected, respectively, at the
start of the epidemic.
(Hence \code{S0 + I0} should be less than or equal to 1.)}
}
\value{
A numeric vector listing values in the interval [0,1] for the expected
epidemic final size.
}
\description{
Computes the proportion of a population expected to be infected over
the course of an epidemic, as a function of the basic reproduction number.
}
\details{
\code{R0}, \code{S0}, and \code{I0} are recycled to length
\code{max(lengths(list(R0, S0, I0)))}.

At least one of \code{S0} and \code{I0} must be supplied.
If \code{S0} (\code{I0}) is supplied but not \code{I0} (\code{S0}),
then the latter is assigned the value of one minus the former.
}
\section{Computation}{

The basic reproduction number \code{R0} defines the expected
epidemic final size \code{Z} through the implicit relation

\code{Z = S0 * (1 - exp(-R0 * (Z + I0)))}.

\code{Z} can be expressed as an explicit function of \code{R0} using the
\href{https://en.wikipedia.org/wiki/Lambert_W_function}{Lambert W function}:

\code{Z(R0, S0, I0) = S0 + (1 / R0) * lambertW(-R0 * S0 * exp(-R0 * (S0 + I0)))}.

\code{finalsize} evaluates this function at each supplied
\code{(R0, S0, I0)} triple.
}

\examples{
R0 <- 10^seq(-3, 1, length.out = 151L)

plot(R0, finalsize(R0, S0 = 1, I0 = 0), type = "l", las = 1,
     xlab = "basic reproduction number",
     ylab = "final size")

}
\references{
Ma J, Earn DJD. Generality of the final size formula for an epidemic
of a newly invading infectious disease. Bull Math Biol. 2006;68:679–-702.
}
\seealso{
\code{\link{timescale}}, \code{\link{R0}}
}
