\name{egf.simulate.egf_model}
\alias{egf.simulate.egf_model}
\title{Estimate Generative Model from Simulation}
\description{
Estimates the generative model for a set of simulated incidence
time series by evaluating a suitable call to \code{\link{egf.egf_model}}.
}
\usage{
\method{egf}{simulate.egf_model}(model, \dots)
}
\arguments{
\item{model}{a \code{\link{simulate.egf_model}} object supplying
  simulated incidence time series and details about the generative model.}
\item{\dots}{optional arguments passed to \code{\link{egf.egf_model}}.}
}
\value{
An \code{\link{egf}} object.
}
\examples{
model <- egf_cache("simulate-egf_model-1.rds")
object <- egf_cache("egf-simulate-egf_model-1.rds", egf(model))
(pp <- cbind(actual = coef(model), fitted = coef(object)))
}
