\name{NEWS}
\title{News for \R{} Package \pkg{epigrowthfit}}

\section{Changes in version 0.15.4 (2025-01-24)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Tests in \file{tests/egf_utils.R} wrongly relied on the
      type of the value of a call to \code{seq.int}, causing errors
      under \R{} 4.5.
    }
  }
}

\section{Changes in version 0.15.3 (2024-06-17)}{
  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item \code{egf_cache} is removed.  Serializing objects under
      \code{tools::R_user_dir("epigrowthfit", "cache")} for reuse in
      examples and tests introduced a race condition when running two
      or more checks concurrently, and workarounds were not appealing.
    }
  }
  \subsection{INSTALLATION}{
    \itemize{
      \item \pkg{emdbook} is removed from the \code{Suggests} field of
      \file{DESCRIPTION}.  A full check of the package requires only
      \pkg{TMB} and its dependencies.
    }
  }
}

\section{Changes in version 0.15.2 (2024-04-24)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Use \verb{\donttest} for examples in
      \file{man/profile.egf.Rd}, which remained too long.

      \item Include URL in \file{man/covid19.ontario.Rd} with
      \verb{\verb} instead of \verb{\href} to avoid HTTP 403.
    }
  }
}

\section{Changes in version 0.15.1 (2024-04-24)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Expand \samp{Description} field in \file{DESCRIPTION}.

      \item Add \verb{\value} tags to two \file{man/*.Rd}.

      \item Clean up after \code{par} calls in \code{plot} method.

      \item Do not assign to global environment in \code{simulate}
      methods. % revert after escaping newbies queue
    }
  }
}

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
