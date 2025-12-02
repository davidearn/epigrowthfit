\name{NEWS}
\title{News for \R{} Package \pkg{epigrowthfit}}
\newcommand{\GH}{\href{https://github.com/davidearn/epigrowthfit/issues/#1}{GH##1}}
\section{Changes in version 0.15.5 (2025-12-02)}{
  \subsection{Description}{
    \itemize{
      \item .
    }
  }
  \subsection{Compatibility}{
    \itemize{
      \item .
    }
  }
  \subsection{Installation}{
    \itemize{
      \item .
    }
  }
  \subsection{Documentation}{
    \itemize{
      \item .
    }
  }
  \subsection{Significant User-Visible Changes}{
    \itemize{
      \item .
    }
  }
  \subsection{Deprecated and Defunct}{
    \itemize{
      \item .
    }
  }
  \subsection{New Features}{
    \itemize{
      \item Function \code{egf1} wrapping function \code{egf} and
      handling its simplest usage, when there is exactly one time series
      and exactly one fitting window.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
  \subsection{Internals}{
    \itemize{
      \item .
    }
  }
}
\section{Changes in version 0.15.4 (2025-02-14)}{
  \subsection{Compatibility}{
    \itemize{
      \item Tests in \file{tests/egf_utils.R} wrongly relied on the
      type of the value of a call to \code{seq.int}, causing errors
      under \R{} 4.5.
    }
  }
}
\section{Changes in version 0.15.3 (2024-06-17)}{
  \subsection{Description}{
    \itemize{
      \item \pkg{emdbook} is removed from \code{Suggests}.  A full check
      of the package requires only \pkg{TMB} and its dependencies.
    }
  }
  \subsection{Deprecated and Defunct}{
    \itemize{
      \item \code{egf_cache} is removed.  Serializing objects under
      \code{tools::R_user_dir("epigrowthfit", "cache")} for reuse in
      examples and tests introduced a race condition when running two
      or more checks concurrently, and workarounds were not appealing.
    }
  }
}
\section{Changes in version 0.15.2 (2024-04-24)}{
  \subsection{Documentation}{
    \itemize{
      \item \verb{\donttest} is used for examples in
      \file{man/profile.egf.Rd}, which remained too long.

      \item \verb{\verb} is used instead of \verb{\href} for URL in
      \file{man/covid19.ontario.Rd} to avoid HTTP 403.
    }
  }
}
\section{Changes in version 0.15.1 (2024-04-24)}{
  \subsection{Description}{
    \itemize{
      \item \samp{Description} field is expanded.
    }
  }
  \subsection{Documentation}{
    \itemize{
      \item \verb{\value} tags are added to two \file{man/*.Rd}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item Methods for \code{plot} clean up after calls to \code{par}.

      \item Methods for \code{simulate} no longer assign to the global
      environment.
    }
  }
}
\section{Changes in version 0.15.0 (2024-04-23)}{
  \subsection{Significant User-Visible Changes}{
    \itemize{
      \item For as long as the version number is 0.y.z, details
      of the API are subject to change without formal deprecation.
      Significant user-visible changes will be documented here.
    }
  }
  \subsection{New Features}{
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
    }
  }
}

%% \section{Changes in version x.y.z (yyyy-mm-dd)}{
%%   \subsection{Description}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Compatibility}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Installation}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Documentation}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Significant User-Visible Changes}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Deprecated and Defunct}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{New Features}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Bug Fixes}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%%   \subsection{Internals}{
%%     \itemize{
%%       \item .
%%     }
%%   }
%% }
