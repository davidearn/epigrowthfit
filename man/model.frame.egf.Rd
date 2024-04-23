\name{model.frame.egf}
\alias{model.frame.egf}
\alias{model.frame.egf_no_fit}
\title{Extract Model Frames}
\description{
Extracts from a model object any of several data frames used to specify
the model, including the mixed effects model frames.
}
\usage{
\method{model.frame}{egf}(formula,
            which = c("ts", "windows", "parameters", "extra", "combined"),
            top = egf_top(formula), full = FALSE, \dots)
}
\arguments{
\item{formula}{
  an \code{\link{egf}} object.}
\item{which}{
  a character string controlling what is returned:
  \describe{
  \item{\code{"ts"}}{
    disease incidence time series.}
  \item{\code{"windows"}}{
    fitting window endpoints.}
  \item{\code{"parameters"}}{
    the mixed effects model frame corresponding to \code{top}.}
  \item{\code{"extra"}}{
    variables preserved in \code{formula} due to setting of
    \code{\link{egf}} argument \code{select_windows}.}
  \item{\code{"combined"}}{
    the result of concatenating (in the sense of \code{cbind})
    all mixed effects model frames and the data frame corresponding
    to \code{"extra"}, then deleting any duplicated variables.}
  }}
\item{top}{
  a character string specifying a top level nonlinear model parameter,
  for \code{which = "parameters"}.}
\item{full}{
  a logical, for \code{which = "ts"}.  If \code{TRUE}, then complete
  time series are returned.  Otherwise, only observations belonging
  to fitting windows are returned.}
\item{\dots}{
  unused optional arguments.}
}
\value{
A data frame.
}
\seealso{
The generic function \code{\link{model.frame}}.
}
