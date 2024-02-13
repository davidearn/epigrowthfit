\name{fitted.egf}
\alias{fitted.egf}
\alias{fitted.egf_no_fit}
\title{Extract fitted values}
\usage{
\method{fitted}{egf}(
  object,
  top = egf_top(object),
  link = TRUE,
  se = FALSE,
  subset = NULL,
  select = NULL,
  ...
)

\method{fitted}{egf_no_fit}(
  object,
  top = egf_top(object),
  link = TRUE,
  se = FALSE,
  subset = NULL,
  select = NULL,
  ...
)
}
\arguments{
\item{object}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{top}{A subset of \code{\link{egf_top}(object, link = TRUE)}
naming top level nonlinear model parameters for which fitted
values should be retrieved.}

\item{link}{A logical flag.
If \code{FALSE}, then fitted values are inverse link-transformed.}

\item{se}{A logical flag. If \code{se = TRUE} and \code{link = TRUE},
then approximate delta method standard errors on fitted values
are reported.
Standard errors are required for subsequent use
of \code{\link{confint.fitted.egf}}.
Setting \code{se = TRUE} and \code{link = FALSE} is an error,
as standard errors are not available for inverse link-transformed
fitted values.
Setting \code{se = TRUE} when \code{object} inherits from class
\code{"egf_no_fit"} is likewise an error.}

\item{subset, select}{Index vectors for the rows and columns of
\code{mf = \link[=model.frame.egf]{model.frame}(object, "combined")}
or language objects evaluating to such vectors.
\code{subset} indicates fitting windows for which confidence intervals
should be computed; the default (\code{NULL}) indicates all of them.
\code{select} indicates variables that should be appended to the result;
the default (\code{NULL}) is to append nothing.
Evaluation of language objects follows the implementation of
\code{\link{subset}(mf, subset, select)}.}

\item{...}{Unused optional arguments.}
}
\value{
A data frame inheriting from class \code{"fitted.egf"}, with variables:
\item{top}{
  Top level nonlinear model parameter,
  from \code{\link{egf_top}(object, link = link)}.
}
\item{ts}{
  Time series, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)$ts)}.
}
\item{window}{
  Fitting window, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)$window)}.
}
\item{estimate}{
  Fitted value of parameter \code{top} in fitting window \code{window}.
  If \code{object} inherits from class \code{"egf"},
  then this is conditioned on the estimated mixed effects model.
  If \code{object} inherits from class \code{"egf_no_fit"},
  then this is conditioned on the \emph{initial} values
  of mixed effects model parameters.
}
\item{se}{
  Approximate delta method standard error on \code{estimate}.
  Absent except for calls matching \code{fitted(link = TRUE, se = TRUE)}.
}
}
\description{
Retrieves fitted values of top level nonlinear model parameters.
The fitted value of a given parameter for a given fitting window
is obtained by adding
(i) the population fitted value computed as a linear combination
of fixed effects coefficients and
(ii) all applicable random effects, with random effects set equal
to their conditional modes.
}
\examples{
object <- egf_cache("egf-1.rds")
zz <- egf_cache("fitted-egf-1.rds", {
    fitted(object, se = TRUE)
})
str(zz)

}
\seealso{
\code{\link{confint.fitted.egf}}

Other extractors: 
\code{\link{coef.egf}()},
\code{\link{df.residual.egf}()},
\code{\link{fixef.egf}()},
\code{\link{formula.egf}()},
\code{\link{getCall.egf}()},
\code{\link{logLik.egf}()},
\code{\link{model.frame.egf}()},
\code{\link{model.matrix.egf}()},
\code{\link{nobs.egf}()},
\code{\link{ranef.egf}()},
\code{\link{terms.egf}()},
\code{\link{vcov.egf}()}
}
\concept{extractors}
