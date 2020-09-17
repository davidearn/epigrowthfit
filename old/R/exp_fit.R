#' exponential fit
#' @param t time vector
#' @param X vector of cases/deaths/etc.
#' @param theta0 A named vector of initial parameter values
#' @param model the model to fit to
#' @param loglik the log-likelihood function
#' @param confint (logical) generate confidence interval?
#' @param lower lower bounds (see \code{link{bbmle::mle2}})
#' @param upper upper bounds (see \code{link{bbmle::mle2}})
#' @param optimizer name of optimization function, e.g. \sQuote{optim} or \sQuote{nlminb} (see \code{\link{mle2}})
#' @param optimfun user-specified optimization function (see \code{\link{mle2}})
#' @param method name of optimization method (ignored unless \code{optimizer} is "optim")
#' @param fallback_stderr_fac quasi-standard error (for profiling), as a proporition of estimated r value, when Hessian is non-positive-definite
#' @param drop_mle2_call to save space, drop call slot from mle2 object
#' 
#' @importFrom stats dpois dnbinom optim as.formula
#' @importFrom stats4 confint profile coef
#' @importFrom bbmle mle2 parnames<- AICc logLik vcov
#' @importFrom utils tail
#' @export
#' @examples
#' dd <- subset(london_bills,outbreak.year==1665)
#' ff <- exp_fit(t=dd$time[15:40], X=dd$plague.deaths[15:40],
#'         theta0=c(r=3,x0=-10,K=10,ll.k=0),
#'         model=get_model("logistic"),
#'         loglik=get_loglik("nbinom"))
#' ffb <- exp_fit(t=dd$time[15:40], X=dd$plague.deaths[15:40],
#'         theta0=c(r=3,x0=-10,K=10,ll.k=0,b=0),
#'         model=get_model("logistic",add_const=TRUE),
#'         loglik=get_loglik("nbinom"),
#'         optimizer="nlminb",
#'         fixed = list(b=-1))
exp_fit <- function(t, X, theta0,
                    model,
                    loglik = get_loglik("nbinom"),
                    confint=TRUE,
                    fixed=NULL,
                    lower=NULL,
                    upper=NULL,
                    optimizer="user",
                    optimfun=NULL,
                    method="L-BFGS-B",
                    optCtrl=list(),
                    extra_optim=FALSE,
                    fallback_stderr_fac = 0.1,
                    drop_mle2_call = TRUE,
                    verbose=FALSE) {
  # the model is for cumulative cases/deaths,
  # thus, the mean is a different model(t1)-model(t2)
  t0 = t[1]
  t1 <- t - t0
  t2 = c(t[-1], tail(t, n=1)+t[2]-t[1]) - t0

  ## save original model for back-transformation later
  model0 <- model

  ## set up model to compute incidence (difference between
  ##   cumulative values computed at time 1 (t1) and time 2 (t2)
  # rename t to t1
  input <- model@input
  model@input <- c(input, "t1", "t2")
  tr <- call("~", input, as.name("t1"))
  model1 <- TransformedModel(model, list(as.formula(tr)))
  # rename t to t2
  tr <- call("~", input, as.name("t2"))
  model2 <- TransformedModel(model, list(as.formula(tr)))
  f <- call("~", as.name(model@output),
            call("-", model2@expr[[1]], model1@expr[[1]]))
  model <- Model(model@name, as.formula(f), input=c("t1", "t2"), par=model@par)
  loglik@input <- c(loglik@input, model@input)
  loglik@par <- names(theta0)
  tr <- call("~", as.name(loglik@input[[1]]), model@expr[[1]])
  loglik <- TransformedModel(loglik, list(as.formula(tr)))

  ## inherits distribution (loglik) from environment ...
  likfun <- function(p) {
    -sum(Eval(loglik, par=p, X=X, t1=t1, t2=t2))
  }
  parnames(likfun) <- names(theta0)

  # the gradient of the likfun
  gr <- function(p) {
    -colSums(grad(loglik, par=p, X=X, t1=t1, t2=t2))
  }
  parnames(gr) <- names(theta0)

  # compute the hessian
  hes <- function(p) {
    -colSums(hessian(loglik, par=p, X=X, t1=t1, t2=t2))
  }
    parnames(hes) <- names(theta0)

    ## default fitting strategy:
    ## 0. custom 'optimizer' optim_combo is actually a combination
    ##   of BFGS and Nelder-Mead, repeating the cycle up to 25 times
    ##   until convergence is achieved, or successive fits have
    ##   log-likelihood decreasing by less than 1e-4
    ## 1. use optim_combo in mle2
    ## 2. after fitting, use optimfun again (why???)
    ## 3. find confidence intervals
    ##    if confint returns something other than 'profile.mle2',
    ##     we need to refit the model (return to #1) with the better
    ##     solution as a starting point, because the parameter
    ##     names in the returned better fit can be wrong (FIXME: why??)
    tries <- 0 
    while (TRUE) {
        tries <- tries+1
        if (optimizer=="user" && is.null(optimfun)) {
            optimfun <- optim_combo
        }
        if (is.null(lower)) lower <- -Inf
        if (is.null(upper)) upper <- Inf
        if (!is.null(fixed)) {
            theta0 <- theta0[!names(theta0) %in% names(fixed)]
        }
        mle2_args <- list(likfun,
                          start=theta0,
                          gr=gr,
                          fixed=fixed,
                          optimizer=optimizer,
                          optimfun=optimfun,
                          control=optCtrl,
                          method=method)
        ## only use strict with optim_combo ...
        if (optimizer=="user") {
            mle2_args <- c(mle2_args,list(strict=TRUE))
        } else if (optimizer=="nlminb") {
            ## if we use lower, upper in optim it may *silently* switch
            ##  to L-BFGS-B?
            mle2_args <- c(mle2_args,list(lower=lower,upper=upper))
        }
        if (verbose) {
            cat(sprintf("fitting mle2 (try #%d)\n",tries))
            if (tries>1) {
                print(theta0)
            }
        }
        fit <- do.call(mle2,mle2_args)
        fit@call$strict <- NULL ## FALSE
        p = coef(fit)

        ## FIXME: how important is this part?
        if (extra_optim && (!is.null(optimfun)) && is.null(fixed)) {
            out <- optimfun(p, likfun, gr, hessian = TRUE, strict = TRUE)
            fit@coef <- fit@fullcoef <- out$par
            ## FIXME: figure out why out$hessian is sometimes NULL??
            try(fit@vcov <- solve(out$hessian),silent=TRUE)
        }
        if (confint) {
            if (verbose) {
                cat(sprintf("computing confint\n"))
            }
            if (is.valid(vcov(fit))) {
                if (verbose) {
                cat(sprintf("profile (easy)\n"))
                }
                prof <- profile(fit, alpha=0.05, which='r',continuation="naive")
            } else {
                if (verbose) {
                    cat(sprintf("profile (hard)\n"))
                }
                prof <- profile(fit, alpha=0.05, which='r',
                                   std.err=p['r']*fallback_stderr_fac,
                                   maxsteps=1000,
                                   continuation="naive")
            }
            if (class(prof) == 'profile.mle2') {
                conf = confint(prof)
                break
            }
            ## if we reach here, confint must have found a 
            ## better solution. we need to repeat
            theta0 = coef(prof)
            names(theta0) = names(p)
        } else break
    }
    if (confint) {
        p.upper <- p.lower <- p
        p.lower[["r"]] <- conf[[1]] 
        p.upper[["r"]] <- conf[[2]]
        result <- c(
            growth.rate = transformPar(model0, p)$r,
            lower = transformPar(model0, p.lower)$r,
            upper = transformPar(model0, p.upper)$r
        )
    } else {
        prof <- NA
        result <- c(growth.rate = transformPar(model0, p)$r,
                    lower=NA, upper=NA)
    }
    if (drop_mle2_call) {
        fit@call <- fit@call.orig <- quote(call_deleted)
    }
    L <- list(result = result, fit = fit, prof = prof, 
              mean = function(p, t) {
                      t1 <- t - t0
	              t2 <- c(t[-1], tail(t, n=1)+t[2]-t[1]) - t0
	              Eval(model, par=p, t1=t1, t2=t2)
	            })
    ## number of re-profiling attempts
    attr(L,"tries") <- tries
    return(L)
}

# a custom optimizing function
#
# make a BFGS try, if not converged to a minimum, try a Nelder-Mead after it.
# Its interface is exactly identical to optim.
#' @importFrom stats na.omit cor
#' @importFrom methods is
optim_combo <- function(par, fn, gr = NULL, ...,
                     method = NULL, lower = -Inf, upper = Inf,
                     control = NULL, hessian = FALSE, strict = TRUE,
                     max.tries=25) {
  L <- Inf
  tries <- 0
  while (tries<max.tries) {
      tries <- tries+1
      out <- optim(par, fn=fn, gr=gr, ...,
                 method='BFGS', lower=-Inf, upper=Inf,
                 control = list(maxit=1e3), hessian = TRUE)
    ## if we made no progress, break
    if (!strict && L - out$value < 1e-4) break
    L <- out$value
    if (out$convergence == 1) {
        warning("optim-not-enough-evaluations")
        ## optim() thinks convergence not reached, keep trying
        ## FIXME: how will it come out differently next time??
        ## next
    }
    # check if the convergence has been reached 
    # and if the vcov is positive definite
    if (out$convergence == 0 && is.valid(out$hessian) &&
        is.valid(solve(out$hessian))) {
        break
    }
    # if we reach here, there is something wrong with the convergence
    # try a Nelder-Mead step
    out <- optim(out$par, fn, method='Nelder-Mead',
                 control = list(maxit=1e5))
    if (isTRUE(all.equal(par,out$par,tolerance=1e-6))) {
      warning("some convergence issue, but optimfun made no progress in finding a minimum. Stop the iteration.")
      break
    }
    par <- out$par
  }
  out
}


##' repeat nlminb until converged
##' @param par parameter
##' @param fn function
##' @param gr gradient
##' @param method stub
##' @param lower lower bound
##' @param upper upper bound
##' @param control control arguments
##' @param hessian stub
##' @param strict ?
##' @param max.tries max number of attempts
##' @export
optim_repeat_nlminb <- function(par, fn, gr = NULL, ...,
                         method = NULL, lower = -Inf, upper = Inf,
                         control = list(), hessian = FALSE, strict = TRUE,
                        max.tries=25) {
    ## for now, hard-code nlminb
    L <- Inf
    tries <- 0
    while (tries<max.tries) {
        tries <- tries+1
        out <- nlminb(start=par, objective=fn, gradient=gr,
                      ...,
                      control = control,
                      lower=lower, upper=upper)
        hessian <- numDeriv::hessian(fn, out$par)
        ## if we made no progress, break
        par <- out$par
        value <- out$objective
        convergence <- out$convergence
        message <- out$message
        if (!strict && L - out$objective < 1e-4) break
        L <- out$objective
        if (out$convergence == 0 && is.valid(hessian) &&
            is.valid(solve(hessian))) {
            break
        }

    }
    return(list(par=par,value=value,convergence=convergence,
                message=message))
}

#' check whether a variance-covariance matrix is valid
#' @param V the vcov matrix
#' @param tol tolerance for positive definiteness
#' @return boolean
is.valid <- function(V,tol=1e-10) {
  if (any(!is.finite(V))) return(FALSE)
  lambda <- eigen(V, symmetric=TRUE, only.values = TRUE)$values
  # all eigenvalues of V must be real and positive (positive definite)
  # all variances must be positive
  all(diag(V) > 0) && all(!is.complex(lambda)) && all(lambda > tol)
}
