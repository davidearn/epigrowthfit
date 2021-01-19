% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.egf}
\alias{predict.egf}
\title{Compute predicted incidence}
\usage{
\method{predict}{egf}(object, subset = NULL, time = NULL, se = FALSE, ...)
}
\arguments{
\item{object}{An \code{"egf"} object returned by \code{\link[=egf]{egf()}}.}

\item{subset}{A named list of atomic vectors of length 1 specifying exactly
one level for each factor in \code{object$frame} (and thus a unique
fitting window). Use the default (\code{NULL}) if and only if
\code{object$frame} has no factors.}

\item{time}{A numeric vector listing increasing time points in days since
the start of the fitting window specified by \code{subset}.}

\item{se}{A logical scalar. If \code{TRUE}, then standard errors on predicted
values (log scale) are also reported.}

\item{...}{Unused optional arguments.}
}
\value{
A list inheriting from class \code{"egf_predict"}, containing two data
frames, \code{log_cum_inc} and \code{log_int_inc}, each with numeric variables
\code{time}, \code{estimate}, and \code{se} (if \code{se = TRUE}).

\code{log_cum_inc$estimate[i]} stores log cumulative incidence at
\code{time[i]}. For \code{i > 1}, \code{log_int_inc$estimate[i]} stores log
incidence during the interval from \code{time[i-1]} to \code{time[i]}.
}
\description{
Computes predicted values of log cumulative incidnece
and log interval incidence given user-specified time points.
}
\details{
Elements of \code{subset} (if non-\code{NULL}) must have the form
\code{factor_name = level_name}, where \code{factor_name} is
the name of a factor in \code{object$frame} and \code{level_name}
is an element of \code{levels(object$frame$factor_name)}.
}
\seealso{
\code{\link[=confint.egf_predict]{confint.egf_predict()}}
}