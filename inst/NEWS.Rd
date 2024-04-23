\name{NEWS}
\title{News for \R{} Package \pkg{epigrowthfit}}

\section{Changes in version 0.15.0 (2024-04-23)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item A \dQuote{traditional} model estimating function \code{egf},
      returning model objects of class \code{egf}.

      \item Methods for generic functions: \code{coef}, \code{confint},
      \code{df.residual}, \code{extractAIC}, \code{fitted},
      \code{fixef} (from \pkg{nlme}), \code{formula}, \code{getCall},
      \code{logLik}, \code{model.frame}, \code{model.matrix},
      \code{nobs}, \code{plot}, \code{predict}, \code{print},
      \code{profile}, \code{ranef} (from \pkg{nlme}), \code{simulate},
      \code{summary}, \code{terms}, \code{vcov}.

      \item Auxiliary functions of use in epidemiology:
      \code{R0}, \code{finalsize}, \code{timescale}, \code{[dpqr]gi}.

      \item For as long as the version number is 0.y.z, details
      of the API are subject to change without formal deprecation.
      Significant user-visible changes will be documented here.
    }
  }
}
