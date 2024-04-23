\name{egf_model}
\alias{egf_model}
\title{Define a Top Level Nonlinear Model}
\description{
Sets flags defining a top level nonlinear model of epidemic growth
to be estimated by \code{\link{egf}}.
}
\usage{
egf_model(curve = c("logistic", "richards", "exponential",
                    "subexponential", "gompertz"),
          excess = FALSE,
          family = c("nbinom", "pois"),
          day_of_week = FALSE)
}
\arguments{
\item{curve}{
  a character string specifying a model for expected cumulative
  incidence as a function of time.}
\item{excess}{
  a logical flag.  If \code{TRUE}, then a constant baseline mortality
  rate is estimated.}
\item{family}{
  a character string specifying a family of discrete probability
  distributions assigned to observations, which are
  the first order differences of observed cumulative incidence.}
\item{day_of_week}{
  an integer flag.  If positive, then day of week effects are
  estimated as offsets relative to the indicated day of week
  (1=Sunday, 2=Monday, and so on).}
}
\value{
A list inheriting from class \code{egf_model} containing the
validated arguments.
}
\seealso{
\code{\link{simulate.egf_model}}.
}
\examples{
model <- egf_model()
str(model)
}
