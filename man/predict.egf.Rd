\name{predict.egf}
\alias{predict.egf}
\alias{confint.predict.egf}
\title{Compute Predicted Values}
\description{
Computes predicted values of top level nonlinear model parameters.
These are conditional on an estimated nonlinear mixed effects model
and, optionally, new data.
}
\usage{
\method{predict}{egf}(object, newdata = NULL, class = FALSE, se = FALSE, \dots)
}
\arguments{
\item{object}{an \code{\link{egf}} object.}
\item{newdata}{a data frame containing variables to replace those
  in the model frame.  The default is to use the model frame as is,
  and currently that is the only implemented behaviour.}
\item{class}{a logical.  If \code{TRUE}, then the value of the
  \code{\link{predict}} call is a \code{\link{predict.egf}} object,
  not a matrix.}
\item{se}{a logical.  If \code{TRUE} and if \code{class = TRUE},
  then the result is augmented with approximate delta method
  standard errors.}
\item{\dots}{additional arguments passed from or to other methods.}
}
\value{
A numeric matrix containing predicted values.

Users can pass \code{class = TRUE} to obtain an augmented result.
Thus, alternatively:

A data frame inheriting from class \code{predict.egf}, with variables:

\item{top}{top level nonlinear model parameter, from
  \code{\link{egf_top}(object)}.}
\item{ts}{time series, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)[["ts"]])}.}
\item{window}{fitting window, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)[["window"]])}.}
\item{value}{predicted value.}
\item{se}{approximate delta method standard error (only if requested).}
}
