\name{egf.egf_model_simulate}
\alias{egf.egf_model_simulate}
\title{Fit simulated incidence time series}
\usage{
\method{egf}{egf_model_simulate}(model, ...)
}
\arguments{
\item{model}{An \code{"\link[=simulate.egf_model]{egf_model_simulate}"} object
supplying simulated time series and details about the generative model
to be estimated.}

\item{...}{Optional arguments passed to \code{\link{egf.egf_model}}.}
}
\value{
An \code{"\link{egf}"} object.
}
\description{
Estimates the generative model underlying a set of simulated
incidence time series.
}
\examples{
model <- egf_cache("simulate-egf_model-1.rds")
object <- egf_cache("egf-egf_model_simulate-1.rds", egf(model))
pp <- cbind(actual = coef(model), fitted = coef(object))

}
