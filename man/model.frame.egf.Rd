\name{model.frame.egf}
\alias{model.frame.egf}
\alias{model.frame.egf_no_fit}
\title{Extract model frames}
\usage{
\method{model.frame}{egf}(
  formula,
  which = c("ts", "windows", "parameters", "extra", "combined"),
  full = FALSE,
  top = egf_top(formula),
  ...
)

\method{model.frame}{egf_no_fit}(
  formula,
  which = c("ts", "windows", "parameters", "extra", "combined"),
  full = FALSE,
  top = egf_top(formula),
  ...
)
}
\arguments{
\item{formula}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{which}{A character string controlling what is returned:
\code{"ts"}, disease incidence time series;
\code{"windows"}, fitting window endpoints;
\code{"parameters"}, the mixed effects model frame
corresponding to parameter \code{top};
\code{"extra"}, variables preserved in \code{formula}
via \code{\link{egf}} argument \code{select_windows}; or
\code{"combined"}, the result of concatenating
(in the sense of \code{\link{cbind}})
the results of \code{"parameters"} and \code{"extra"},
then deleting any duplicated variables.}

\item{full}{A logical flag.
If \code{TRUE}, then complete time series are returned.
Otherwise, only observations belonging to fitting windows are returned.
Ignored if \code{which != "ts"}.}

\item{top}{A character string specifying a top level nonlinear model parameter.
Ignored if \code{which != "parameters"}.}

\item{...}{Unused optional arguments.}
}
\value{
A data frame.
}
\description{
Extract data frames, including \link[=model.frame]{model frame}s,
from a model object.
}
\details{
Only \code{which = "parameters"} results in a "proper" model frame,
i.e., a data frame with a \code{\link[=terms.object]{terms}} attribute.
}
\seealso{
Other extractors: 
\code{\link{coef.egf}()},
\code{\link{df.residual.egf}()},
\code{\link{fitted.egf}()},
\code{\link{fixef.egf}()},
\code{\link{formula.egf}()},
\code{\link{getCall.egf}()},
\code{\link{logLik.egf}()},
\code{\link{model.matrix.egf}()},
\code{\link{nobs.egf}()},
\code{\link{ranef.egf}()},
\code{\link{terms.egf}()},
\code{\link{vcov.egf}()}
}
\concept{extractors}
