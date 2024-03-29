\name{egf_model}
\alias{egf_model}
\title{Define an epidemic model}
\usage{
egf_model(
  curve = c("logistic", "richards", "exponential", "subexponential", "gompertz"),
  excess = FALSE,
  family = c("nbinom", "pois"),
  day_of_week = FALSE
)
}
\arguments{
\item{curve}{A character string specifying a model for expected cumulative
disease incidence as a function of time.}

\item{excess}{A logical flag. If \code{TRUE}, then a constant baseline mortality rate
is estimated. Set to \code{TRUE} if what is observed is multiple causes
mortality rather than disease mortality or disease incidence.}

\item{family}{A character string specifying a family of discrete probability
distributions assigned to observations of disease incidence.}

\item{day_of_week}{An integer flag. If positive, then day of week effects are estimated
as offsets relative to the indicated day of week
(Sunday if \code{day_of_week = 1}, Monday if \code{day_of_week = 2},
and so on). Logical values are coerced to integer.}
}
\value{
A list inheriting from class \code{"egf_model"} containing the arguments
(after possible matching and coercion).
}
\description{
Sets flags defining the top level nonlinear model of epidemic growth
to be estimated by \code{\link{egf}}.
}
\examples{
model <- egf_model()
str(model)

}
\seealso{
Other behaviour-defining functions: 
\code{\link{egf_control}()},
\code{\link{egf_optimizer}()},
\code{\link{egf_parallel}()}
}
\concept{behaviour-defining functions}
