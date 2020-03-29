#' epidemic growth fit
#'
#' convert a julian date to a date, being careful about truncation (\code{\link{as.Date}} truncates values instead of rounding)
#' @param tvec vector of julian dates (YYYY.fff)
#' 
#' @export
safe_date <- function(tvec) {
    ##  as.Date() **TRUNCATES**, magnifying possible roundoff error:
    ## round first!
    return(as.Date(round(lubridate::date_decimal(tvec))))
}

#' @description
#' A class representing an epidemic growth fit obtained
#' using the methodology of Ma \emph{et al.} (2014)
#' 
#' @slot model the trajectory model (an object of class \code{\link{Model}})
#' @slot loglik the log-likelihood for the observation error distribution
#' @slot time decimal dates at each point in the time series
#' @slot date time as \code{\link{Date}} object
#' @slot deaths deaths in the time series
#' @slot window the fitting window (indices of times)
#' @slot mle2 an \code{\link{mle2}} object
#' @slot growthRate the fitted exponential growth rate
#' @slot mean predicted mean value?
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#' @include models.R  
#' @importFrom methods setClass
#' @rdname epigrowthfit
#' @name epigrowthfit
#' @export
epigrowthfit <- setClass(
  "epigrowthfit",
  slots = c(
    call = "call",
    model = "Model",
    loglik = "Model",
    time = "numeric",
    date = "Date",
    deaths = "numeric",
    window = "numeric",
    first = "numeric",
    peak = "numeric",
    last = "numeric",
    mean = "function",
    mle2 = "mle2",
    growthRate = "numeric",
    profile = "profile.mle2"
  )
)

#' @importFrom stats smooth.spline
get.smooth.max <- function(count,spar=1) {
    ss <- smooth.spline(1:length(count),log(count+1),spar=spar)
    dd <- predict(ss,deriv=1)$y
    ## change in sign of first derivative
    dds <- diff(sign(dd))
    peakvals <- count[dds<0]
    return(peakvals[1])
}

#' @rdname epigrowthfit
#' @importFrom stats lm na.omit
#' @inheritParams epigrowthfit
#' @export
getInits <- function(time,deaths,
                     theta0=NULL,
                     peak=NULL,
                     first=NULL,
                     first_level=NULL,
                     min_window_len = 5,
                     distrib,
                     model,
                     peak_method="max",
                     skip_zero=TRUE,
                     min_r=0.1,
                     debug_window=FALSE,
                     debug_plot=FALSE) {

    peak_extend <- 1
    ## a scaling factor
    N = sum(deaths)
    if (is.null(peak)) {
        ## make an initial guess
        peak = switch(peak_method,
                      max=which.max(deaths),
                      smooth_max=get.smooth.max(deaths))
    }
    ## pick the start of a fitting window, if not specified
    ## choose the last trough before the peak
    if (is.null(first)) {
        if (is.null(first_level)) {
            ## find trough: defined as minimum value before the peak
            ## which.min finds the *first* such value, so which.min(rev(x))
            ## will find the *last* value
            rtrough = which.min(rev(deaths[1:peak]))
        } else {
            below_min <- deaths[1:peak]<(first_level*max(deaths))
            rtrough = if (!any(below_min)) peak else which(rev(below_min))[1]
        }
        first = peak + peak_extend - rtrough
    }
    wlen <- (peak + peak_extend - first) + 1 ## inclusive window length
    ## we start from a nonzero initial condition
    while(skip_zero && deaths[first] == 0) first = first + 1
    if (debug_window) {
        return(list(first=first,peak=peak,wlen=wlen,
                    min_window_len=min_window_len,
                    peak_extend=peak_extend))
    }
    ## try to satisfy minimum window length specification
    ##  (may interact strangely with non-zero IC criterion below?)
    if (wlen<min_window_len) {
        warning("extending window to meet minimum condition")
        first <- max(1,peak+peak_extend-min_window_len)
    }

    ## replace element nm in x with val if missing or NA
    repl <- function(x,nm,val) {
        if ((n <- length(nm))>1) {
            for (i in seq(n)) {
                x <- repl(x,nm[i],val[i])
            }
        } else {
            if (length(nm)>1 || length(x[nm])>1) stop("oops")
            if (!(nm %in% names(x)) || is.na(x[nm])) {
                x[nm] <- val
            }
        }
        return(x)
    }
        
    if (is.null(theta0) || !all(c("r","K","N") %in% names(theta0))) {
        l = floor((peak + peak_extend - first)/2)
        if (l == 1) l = 2
        dd = data.frame(time,deaths)[first:(first+l),]
        ## add a bit to allow zeros in window
        guess = lm(log(deaths+0.1) ~ I(time-time[1]),dd) 
        p = coef(guess)
        r_est <- p[[2]]
        if (is.na(r_est)) stop("getInits: r_est is NA")
        if (r_est<=0) {
            warning("estimated r<=0: setting to min value ",min_r)
            r_est <- min_r
        }
        if (debug_plot) {
            with(dd,plot(time-time[1],log(deaths+0.1)))
            abline(guess)
        }
        theta0 <- repl(theta0,c("r","x0","K"),c(r_est,exp(p[[1]]),N))
    }
    if (is.character(model)) {
        model <- get_model(model)
    }
    model_name <- model@name
    if (model_name=="richards" && !("s" %in% names(theta0))) {
        theta0 <- c(theta0, s = 1.0001)
    }
    ## negative binomial has an extra parameter
    if (distrib=="nbinom" && (!("ll.k" %in% names(theta0)))) {
        theta0 <- c(theta0,ll.k=1)
    }
    return(list(first=first,peak=peak,theta0=theta0))
}

#' epidemic growth fit
#' @rdname epigrowthfit
#' @name epigrowthfit
#' 
#' @description
#' Fitting function for plague data. Uses either Richards (by default, or if an alpha parameter is supplied in \code{theta0}) or logistic model (if alpha not supplied, or if data set contains <= 4 points).
#'
#' \code{getInits} picks an appropriate fitting window (see Details).
#'
#' @details
#' \describe{
#' \item{fitting window}{if \code{peak} is \code{NULL} and \code{first_level} is non-\code{NULL} (the default), \code{peak} is set to the maximum of the data, while \code{first} is set to the first element after the last observation with cases/deaths less than (\code{first_level*max(data)}). Otherwise, if \code{first_level} and \code{first} are both \code{NULL}, \code{first} is set to the first element after the last local minimum before the peak. If either \code{first} or \code{peak} are non-\code{NULL} (integers), they are set to the specified values. The fitting window goes from elements \code{first} to \code{peak}+1 (inclusive) of the data (or to \code{peak} for exponential fits).}
#' \item{initial parameter values}{If initial parameters (\code{theta0}) are not supplied, uses log-linear fitting over a range starting at \code{first} and going halfway to \code{peak+1} to estimate the initial size and growth rates. The initial population size (K) is set to the total number of cases/deaths. The Richards parameter (s) is initially set to 1.0001; the log of the negative binomial size parameter is initially set to 1.}
#' }
#' \itemize{
#' \item{to debug, use \code{trace("initialize",sig="epigrowthfit",browser)}
#' }
#' }
#' @param time time vector
#' @param deaths death vector
#' @param data data frame
#' @param time_var column name of time vector within \code{data}
#' @param deaths_var column name of deaths vector within \code{data}
#' @param theta0 initial parameter guess
#' @param first first time step to use
#' @param first_level threshold to pick start of the fitting window (as a fraction of the peak incidence)
#' @param peak estimated peak of epidemic (fitting window ends at peak
#' for exponential models, peak+1 otherwise; initial parameter guesses
#' are based on log-linear fits to the data from the first to
#' halfway between the first step #' and the peak)
#' @param first_time start time for fitting window
#' @param last_time end time for fitting window
#' @param baseline A "baseline" refers to a cumulative incidence model that incorporates
#' a \code{b*t} term.  If \code{baseline=TRUE} then (if a baseline is being used) 
#' its value (\code{b}) is fixed (i.e., it is a fixed offset), based on
#' the computed mean values in the surrounding years, rather than estimated
#' from the outbreak year data.  Whether or not a baseline is used is specified
#' by the \code{add_const} argument of the \code{\link{get_model}} function.
#' @param model model to use
#' @param distrib observation error distribution to use (FIXME:
#' slot name should match argument?  (distrib != loglik)
#' @param transforms a list of formulas specifying the parameter transformations, both for the model and for the log-likelihood function. Each transform must be a formula specifying a parameter depending on a smooth function of itself, e.g., x0~exp(x0).
#' @param inverses a list of formulas specifying the inverse parameter transformations, both for the model and for the log-likelihood function. Each transform must be a formula specifying a parameter depending on a smooth function of itself, e.g., x0~exp(x0).
#' @return If \code{return_fit} is \code{FALSE}, a list containing elements \code{growth.rate} (point estimate (\code{value} and lower/upper CIs for the growth rate (\code{lower}, \code{upper})) and \code{fit} (predictions: \code{time}, observed \code{data}, and fitted values \code{fit}). If \code{return_fit} is \code{TRUE}, the results of \code{\link{exp_fit}}
#' @examples
#' ## extract Great Plague of London from full plague data frame
#' great.plague <- subset(london_bills,outbreak.year==1665)
#' ## run default fit, which uses Richards model
#' gpfit.richards <- epigrowthfit(data=great.plague)
#' ## run fit using logistic model
#' gpfit.logistic <- epigrowthfit(data=great.plague,model="logistic")
#' ## run fit using exponential model
#' gpfit.exp <- with(great.plague,epigrowthfit(time,plague.deaths,model="exp"))
#' ## print growth rates and confidence intervals for each fit
#' growthRate(gpfit.richards)
#' growthRate(gpfit.logistic)
#' growthRate(gpfit.exp)
#' ## print associated doubling time, R0 and fitting window range for Richards model fit
#' doublingTime(gpfit.richards)
#' R0(gpfit.richards)
#' range(window(gpfit.richards))
#' gpfit.richards
#' ## spew full summary of this fit
#' summary(gpfit.richards)
#' ## plot fit in various ways
#' plot(gpfit.richards)
#' plot(gpfit.richards,log="y")
#' plot(gpfit.richards,log="y",window.only=TRUE)
#' plot(gpfit.richards,cumulative=TRUE)
#' plot(gpfit.richards,cumulative=TRUE,log="y")
#' ## plot all fits on same graph for comparison
#' plot(gpfit.richards,lwd.pred=6,show.legend=FALSE,show.gof=FALSE)
#' plot(gpfit.logistic,add=TRUE,col.pred="green",lwd.pred=1,show.data=FALSE)
#' plot(gpfit.exp,add=TRUE,col.pred="blue",show.data=FALSE)
#' ## add legend showing goodness of fit of each model
#' legend.text <- sprintf("%s (%.3f)",
#'   c("Richards","logistic","exp"),
#'   c(gof(gpfit.richards),gof(gpfit.logistic),gof(gpfit.exp)))
#' legend("topright",bty="n",legend=legend.text,lty=1,lwd=c(6,1,3),col=c("black","green","blue"))
#' title(main="Initial Growth Fits to Great Plague of London in 1665")
#' \dontrun{
#'     ## hard-code fitting window
#'     gpfit.3 <- update(gpfit.logistic, first_time=1665.5, last_time=1666.0)
#'     ## constant-baseline fits
#'     d1375 <- subset(husting_wills,outbreak.year==1375)
#'     gpfit.4 <- epigrowthfit(data=d1375,
#'           model=get_model("logistic",add_const=TRUE),
#'           first=1,peak_extend=4,
#'           deaths_var="nwills",
#'           baseline=TRUE)
#' }
#' if (require("outbreaks")) {
#' #' ## example using outbreaks package data
#' }
#' @importFrom methods getClass initialize new is
#' @export epigrowthfit
setMethod(
  "initialize",
  "epigrowthfit",
  definition = function(.Object, time=NULL, deaths=NULL,
                        data=NULL,
                        time_var="time",
                        deaths_var="plague.deaths",
                        baseline=NULL,
                        fixed=NULL,
                        theta0=NULL,
                        first_time = NULL,
                        first = NULL,
                        first_level = 0.02,
                        min_window_len = 5,
                        peak = NULL,
                        peak_extend = NULL,
                        last = NULL,
                        last_time = NULL,
                        skip_zero = TRUE,
                        link_vals = NULL,
                        model = getOption("epigrowthfit.model","richards"),
                        distrib=getOption("epigrowthfit.distrib","nbinom"),
                        ## would like to pass these to exp_fit via ...,
                        ## but seems there is some magic here ...
                        optimizer = getOption("epigrowthfit.optimizer","user"),
                        optimfun = getOption("epigrowthfit.optimfun",NULL),
                        method = getOption("epigrowthfit.method","L-BFGS-B"),
                        upper = Inf,
                        optCtrl = getOption("epigrowthfit.optCtrl",list()),
                        verbose=FALSE,
                        ...) {
    ## need to go back 3 steps to find "original" call
    .Object@call <- match.call(call=sys.call(sys.parent(3)))
    ## .Object <- callNextMethod(.Object, ...)
    time <- .Object@time <-
        if (!is.null(time)) time else data[[time_var]]
    deaths <- .Object@deaths <-
        if (!is.null(deaths)) deaths else data[[deaths_var]]
    ## date is convenient for plotting results:
    .Object@date <- safe_date(time)

    loglik <- get_loglik(distrib)
    .Object@loglik <- loglik
    
    ## we use either the Richards model or the logistic model, 
    ## the model is selected by the presence of the alpha parameter 
    ## in the initial guess, representing the Richards model.
    ## or, if data set is too small (=4) we use the logistic model
    ## default to richards
    if (is.null(model)) {
        if (is.null(theta0) || "s" %in% names(theta0)) {
            model <- get_model("richards", link_vals=link_vals)
        } else if ("K" %in% names(theta0)) {
            model <- get_model("logistic")
        } else model <- get_model("exp", link_vals=link_vals)
    }
    if (is.character(model)) {
        model <- get_model(model, link_vals=link_vals)
    }
    model_name <- model@name

    if (is.null(theta0)) {
        theta0 <- c(model@inits,loglik@inits)
    } else {
        xnms <- names(theta0)[is.na(theta0)]
        theta0[xnms] <- c(model@inits[xnms],loglik@inits[xnms])
    }

    ## if specified, set first from first_time **BEFORE** calling getInits ...
    if (!is.null(first_time)) {
        first <- which.min(abs(time-first_time))
    }

    ## FIXME: should pass ... or a list of init options
    inits <- getInits(time,deaths,theta0,
                      peak=peak,
                      first=first,
                      first_level=first_level,
                      distrib=distrib,
                      skip_zero=skip_zero,
                      model=model_name)

    .Object@first <- inits$first

    ## recover results ...
    theta0 <- inits$theta0
    peak <- inits$peak
    first <- inits$first
    ## save in slots:
    .Object@peak <- peak

    ## pick the end of the fitting window
    ## if the time series has data after the peak, we pick the fitting window as
    ## first:(1 step after peak)
    ## except that for the exponential model, we stop at the peak

    if (is.null(last)) {
        if (!is.null(last_time)) {
            last <- which.min(abs(x=time-last_time))
        } else {
            if (peak == length(deaths)) {
                pe <- 0
            } else if (is.null(peak_extend)) {
                pe <- if (model_name=="exp") 0 else 1
            } else pe <- peak_extend
            last <- peak + pe
        }
    }
    .Object@last <- last
    ## check for consistency with data length
    len <- last - first + 1 - length(loglik@par)
    model_has_const <- grepl("b ?\\* ?t",deparse(model@orig_expr))
    if (len == 4 && model_name == " richards") {
        warning("only 4 data points available, using the logistic model instead")
        model <- get_model("logistic", add_const=model_has_const,
                           link_vals=link_vals)
    } else if (len == 3) {
        warning("only 4 data points available, using the exp model instead")
        model <- get_model("exp", add_const=model_has_const,
                           link_vals=link_vals)
    } else if (len < 3)
        stop ("not enough data to fit.", call.=FALSE)
    ## window holds the fitting window
    .Object@window <- window <- first : last

        ## rescale x0 parameter to population scale (FIXME:: document this ...
        ##   and/or make it part of the parameter transformation machinery)
  	if (model_name %in% c("richards","logistic")) {
  	  theta0[["x0"]] <- theta0[["x0"]] / theta0[["K"]]
  	}
  	.Object@model <- model
  	# remove extra parameters
  	theta0 = theta0[names(theta0) %in% c(model@par, loglik@par)]

        ## "inverse" means transform to the link scale ...
  	theta0 <- unlist(transformPar(model, theta0, inverse=TRUE))
  	theta0 <- unlist(transformPar(loglik, theta0, inverse=TRUE))

        if (!is.null(baseline) && !identical(baseline,FALSE)) {
            if (isTRUE(baseline)) {
                baseline <- extract_baseline(data)
            }
            stopifnot(is.numeric(baseline))
            fixed <- list(b=log(unname(baseline)))
        }
    
  	# fit
  	res <- exp_fit(t=time[window]-time[first], 
                       X=deaths[window], 
                       theta0=theta0,
                       fixed=fixed,
                       model=model,
                       loglik=loglik,
                       optimizer=optimizer,
                       optimfun=optimfun,
                       method=method,
                       upper=upper,
                       optCtrl=optCtrl,
                       verbose=verbose,
                       ...)

    .Object@mle2 <- res$fit
    .Object@growthRate <- c(value = res$result[["growth.rate"]],
                            lower = res$result[["lower"]],
                            upper = res$result[["upper"]])
      .Object@mean <- res$mean
      .Object@profile <- res$prof
      coefWarn(.Object)
      ## add data attributes to the fit object:
      data.attr.names <- setdiff(names(attributes(data)),c("names","row.names","class"))
      for (a in data.attr.names) {
          attr(.Object,a) <- attr(data,a)
      }
      ## store starting values
      attr(.Object,"theta0") <- theta0
      .Object
  }
)

#' returns R0 and its 95\% confidence interval
#' @description uses the exponential growth rate from the expgrowthrate
#' object that is passed, and the serial interval distribution available
#' in the package
#' @param obj a fitted model object
#' @docType methods
#' @return a vector
#' @importFrom methods setMethod
#' @export
#' @references
#' \insertRef{WallLips07}{epigrowthfit}
setGeneric(
  "R0",
  def=function(obj, ...) {
    standardGeneric("R0")
  }
)

#' returns R0 and its 95\% confidence interval
#' @param generation.interval the generation interval distribution
#' @docType methods
#' @return a vector
#' @param obj a fitted model object
#' @importFrom methods setMethod
#' @export
#' @references
#' \insertRef{WallLips07}{epigrowthfit}
setMethod(
  "R0",
  "epigrowthfit",
  definition = function(obj, generation.interval=generation.interval.plague()) {
    R0(obj@growthRate, generation.interval)
  }
)

#' returns growth rate and its 95\% confidence interval
#' @docType methods
#' @return a vector
#' @importFrom methods setMethod
#' @export
setGeneric(
    "growthRate",
    def=function(obj, ...) {
    standardGeneric("growthRate")
}
)

#' returns the growth rate and its 95\% confidence interval
#' @docType methods
#' @return a vector
#' @importFrom methods setMethod
#' @export
setMethod(
    "growthRate",
    signature="epigrowthfit",
    definition = function(obj) {
    return(obj@growthRate)
}
)

#' returns the goodness of fit of epigrowthfit object
#'
#' @description
#' goodness of fit is defined as \code{cor(predicted,observed)^2}
#' @docType methods
#' @return a vector
#' @importFrom methods setMethod
#' @importFrom stats cor
#' @export
setGeneric(
  "gof",
  def=function(obj, ...) {
    standardGeneric("gof")
  }
)

#' returns the goodness of fit of epigrowthfit object
#' @docType methods
#' @return a vector
#' @importFrom methods setMethod
#' @export
setMethod(
  "gof",
  "epigrowthfit",
  definition = function(obj,agg=NULL,debug=FALSE) {
    pred <- predict(obj)$deaths
    obs <- obj@deaths[obj@window]
    if (!is.null(agg)) {
        dd <- data.frame(time=obj@time[obj@window],
                         date=obj@date[obj@window],
                         nwills=obj@deaths[obj@window])
        dd_agg <- with(dd,aggsum(deaths=nwills,time=time,period=agg))
        perdiff <- abs(diff(dd_agg$time)[1]-diff(dd$time)[1])
        if (perdiff<1e-4) {
            warning("small agg difference: skipping aggregation steps")
        } else {
            dd_agg$pred <- predict(obj,dd_agg$time)$deaths
            ## predict goes a little bit beyond the aggregated data ... ?
            dd_agg <- na.omit(dd_agg)
            obs <- dd_agg$deaths
            pred <- dd_agg$pred
            if (debug) {
                plot(dd_agg$time,obs)
                lines(dd_agg$time,pred,type="b",pch=4)
                abline(v=range(dd$time),lty=2)
            }
        }
    }
    ## RMSE
    nullvar <- sum((mean(obs)-obs)^2)
    predvar <- sum((pred-obs)^2)
    gofval <- (nullvar-predvar)/(nullvar)
    return(gofval)
  }
)

#' returns the doubling time and its 95\% confidence interval
#' @return a vector
#' @docType methods
#' @importFrom methods setMethod
#' @export
setGeneric(
  "doublingTime",
  def=function(obj) {
    standardGeneric("doublingTime")
  }
)

#' returns the doubling time and its 95\% confidence interval
#' @return a vector
#' @docType methods
#' @importFrom methods setMethod
#' @export
setMethod(
  "doublingTime",
  "epigrowthfit",
  definition = function(obj) {
    v = c(value = obj@growthRate[["value"]],
          lower = obj@growthRate[["upper"]],
          upper = obj@growthRate[["lower"]])
    doublingTime(v)
  }
)


#' returns the fitted values
#' @return a dataframe object with columns time and deaths
#' @importFrom methods setMethod
#' @importFrom stats fitted quantile
#' @importFrom MASS mvrnorm
#' @export
fitted.epigrowthfit <- function(object, CI=TRUE, ...) {
  ## FIXME: how is this different from predict method?
  ## MVN sampling is obsolete/overlaps with predict?  
  if (length(list(...))>0) warning("extra arguments ignored")
  time <- object@time[object@window]
  deaths <- object@mean(coef(object@mle2), time - time[1])
  if (CI) {
      ## sort out fixed/non-fixed parameters
      fc <- coef(object@mle2,exclude.fixed=FALSE)
      fw <- coef(object@mle2,exclude.fixed=TRUE)
      fparind <- which(!names(fc) %in% names(fw))
      fparnms <-  setdiff(names(fc),names(fw))
      runs <- 500
      samples <- matrix(NA, ncol=runs, nrow=length(time))
      V <- vcov(object@mle2)
      if (!is.valid(V,tol=0) && ("ll.k" %in% names(fw))) {
          ## try harder
          ll.k_pos <- which(names(fw)=="ll.k")
          V <- solve(object@mle2@details$hessian[-ll.k_pos,-ll.k_pos],tol=0)
          fw <- fw[-ll.k_pos]
      }
      if (is.valid(V,tol=0)) {
        p <- mvrnorm(runs, mu=fw, Sigma=V)
        if (length(fparnms)>0) {
            p <- cbind(p,fc[fparnms])
            colnames(p)[ncol(p)] <- fparnms
        }
        t <- object@time[object@window]
        for (i in 1:runs) {
            samples[, i ] <- object@mean(p[i,], t-t[1])
        }
      ## ?? ## ind <- which(samples[1,]<0)
      CIs <- apply(samples, 1, function(x)
        quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
    } else CIs=matrix(NA, nrow=2, ncol=length(time))
  }
  d <- data.frame(time=time, deaths=deaths)
  if (CI) cbind(d, lower = CIs[1,], upper = CIs[2,]) else d 
}

#' returns the fitting window
#' @return a vector
#' @docType methods
#' @importFrom methods setMethod
#' @importFrom stats window
#' @param x an \code{epigrowthfit} object
#' @param ... extra arguments (ignored)
#' @export
window.epigrowthfit <- function(x, ...) {
  x@time[x@window]
}

coefWarn <- function(object,tol=10) {
    upar <- coef(object,"unconstrained")
    npar <- coef(object,"natural")
    upar <- upar[upar!=npar]  ## only check non-identity links
    bigcoefs <- abs(upar)>tol
    if (any(bigcoefs)) {
        warning("some extreme parameters estimated: ",
                paste(names(upar)[bigcoefs],collapse=", "))
    }
    invisible(NULL)
}

setMethod("show",
          "epigrowthfit",
          definition = function(object) {
    coefWarn(object)
              cat("EPIDEMIC GROWTH FIT:\n")
              d <- object@date
              cat(sprintf("\t   full date range:\t%s - %s\t[%d points]\n",min(d),max(d),length(d)))
              t <- object@time
              cat(sprintf("\tdecimal date range:\t%.5f - %.5f\n",
                          min(t),max(t)))
              w <- object@window
              cat(sprintf("\tfitted date window:\t%s - %s\t[%d points]\n",min(d[w]),max(d[w]),length(d[w])))
              cat(sprintf("\t     decimal dates:\t%.5f - %.5f\n",
                          min(t[w]),max(t[w])))
              cat(sprintf("\tfitted index window:\t%d:%d\n",
                          min(w),max(w)))
              cat(sprintf("\tfitted peak index:\t%d\n",
                          object@peak))
              cat("\ttrajectory model: ",object@model@name,": ",
                  as.character(object@model@expr),"\n")
              cat("\tdistribution: ",object@loglik@name,"\n")
              cat("\tparameter estimates:",
                  sprintf(" %s = %g", names(coef(object)), coef(object)),
                  "\n")
              cat("\tgoodness of fit [cor(pred,obs)^2]: ", gof(object), "\n")
          })

## return more details
#' @importFrom methods show
#' @export
setClass("summary.epigrowthfit",
         representation(fulltime="numeric",
                        ntime="integer",
                        wintime="numeric",
                        winrange="integer",
                        modelname="character",
                        distribname="character",
                        growthrate="numeric",
                        doublingtime="numeric",
                        R0="numeric",
                        finalsize="numeric",
                        coefficients="numeric",
                        gof="numeric",
                        conv="numeric",
                        convmsg="character")
         )

#' @importFrom stats setNames
#' @export
setMethod("summary",
          "epigrowthfit",
          definition=function(object,gof_agg=NULL,
                              plague_R0=FALSE) {
              rr <- function(x) setNames(range(x),c("min","max"))
              msg <- object@mle2@details$message
              if (is.null(msg)) msg <- ""
    return(new("summary.epigrowthfit",
               fulltime=rr(object@time),
               ntime=length(object@time),
               wintime=rr(object@time[object@window]),
               winrange=rr(object@window),
               modelname=object@model@name,
               distribname=object@loglik@name,
               growthrate=growthRate(object),
               doublingtime=doublingTime(object),
               R0=if (plague_R0) R0(object) else NA_real_,
               finalsize=if (plague_R0) finalsize(R0(object)) else NA_real_,
               coefficients=coef(object),
               gof=gof(object,gof_agg),
               conv=object@mle2@details$conv,
               convmsg=msg))
})

#' @importFrom methods slotNames slot
#' @export
setAs("summary.epigrowthfit", "data.frame",
      function(from) {
    fields <- slotNames(from)
    v <- setNames(lapply(fields,slot,object=from),fields)
    numfields <- sapply(v,is.numeric)
    do.call(data.frame,c(as.list(unlist(v[numfields])),
                         as.list(v[!numfields]),
                         list(stringsAsFactors=FALSE)))
})

#' @export
setMethod("show",
          "summary.epigrowthfit",
          definition=function(object) {
    cat("EPIDEMIC GROWTH FIT SUMMARY:\n")
    cat(sprintf("Full time range:\t%.2f - %.2f\t[%d points]\n",
                object@fulltime[1],object@fulltime[2],object@ntime))
    cat(sprintf("Fitted time window:\t%.2f - %.2f\t[%d points]\n",
                object@wintime[1],object@wintime[2],diff(object@winrange)+1))
    cat(sprintf("Fitted index window:\t%d:%d\n",
                object@winrange[1],object@winrange[2]))
    cat("Trajectory model: ",object@modelname,"\n")
    cat("Distribution: ",object@distribname,"\n")
    cat("Growth rate:\n");  show(object@growthrate)
    cat("Doubling time:\n");  show(object@doublingtime)
    cat("Basic reproduction number R0:\n"); show(object@R0)
    cat("Final size of epidemic:\n"); show(object@finalsize)
    cat("Parameter estimates:\n"); show(object@coefficients)
    cat("Goodness of fit:\n"); show(object@gof)
    return(invisible(object))
})

#' plot method for epigrowthfit objects
#'
#' @rdname plot.epigrowthfit
#' @description
#' simple plotting method to show observed data and fitted model
#' @param x an \code{\link{epigrowthfit}} object
#' @param window.only logical flag: if \code{TRUE} then restrict plot to fitting window
#' @param cumulative logical flag: if \code{TRUE} then plot cumulative data
#' @param add logical flag: if \code{TRUE} then add to existing plot
#' @param show.data logical flag: if \code{TRUE} then display the observed data
#' @param show.pred logical flag: if \code{TRUE} then display the fitted model
#' @param show.legend logical flag: if \code{TRUE} then display trajectory model and noise model in a legend
#' @param show.gof logical flag: if \code{TRUE} then display goodness of fit in a legend
#' @param show.totals logical flag: if \code{TRUE} then display total deaths in all and in the fitting window in a legend
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param log log scale ("", "x", "y", or "xy", as in \code{\link{plot.default}})
#' @param col colour of lines and points for observed data
#' @param col.bg fill colour for observed data points
#' @param col.pred colour of fitted model
#' @param lwd.pred line width of fitted model
#' @param col.year colour of year labels on x axis (in multi-year plots)
#' @param confint include confidence interval ribbon in plot?
#' @param col.ci colour of confidence interval ribbon
#' @param ... further arguments passed to \code{\link{plot}} and \code{\link{points}}
#' @return nothing at the moment
#' @docType methods
#' @importFrom methods setMethod
#' @importFrom grDevices adjustcolor
#' @importFrom graphics lines matplot par
#' @export
setMethod("plot",
          ## plot signature *requires* 'x', 'y' as argument names
          signature(x="epigrowthfit",y="missing"),
          definition = function(x,
                                window.only=FALSE, cumulative=FALSE,
                                show.data=TRUE, show.pred=TRUE, show.legend=TRUE,
                                show.gof=TRUE, show.totals=TRUE,
                                xlab="",ylab=NULL,
                                col="black",
                                col.bg="red",
                                col.pred="black",lwd.pred=3,
                                col.year="blue",
                                confint=FALSE,
                                col.ci=adjustcolor("black",alpha.f=0.3),
                                log="",
                                add=FALSE, ...) {
              w <- if (window.only) x@window else seq_along(x@time)
              time <- x@time[w]
              date <- x@date[w]
              y <- if (cumulative) cumsum(x@deaths[w]) else x@deaths[w]
              total.deaths <- sum(x@deaths[w])
              ## set up plot
              if (!add) {
                  if (is.null(ylab)) {
                      aggregation <- if (is.null(aa <- attr(x,"aggregation"))) "" else aa
                      type <- if (is.null(tt <- attr(x,"type"))) "Deaths" else tt
                      place <- if (is.null(pp <- attr(x,"place"))) "???" else pp
                      ylab <- sprintf("%s %s in %s",aggregation,type,place)
                  }
                  plot(date,y,type="n",las=1,xlab=xlab,ylab=ylab,
                       xaxt="n",log=log,...)
                  df <- data.frame(year=trunc(time),time=time,date=date)
                  fancy.axis.Date(df, col.year=col.year)
              }
              ## indicate the fitting window unless only the window is shown
              if (!window.only)
                abline(v=x@date[range(x@window)],
                       col="grey",lwd=3,lty="dotted")
              ## add the data
              if (show.data) points(date,y,col=col,pch=21,bg=col.bg,type="o",...)
              pred <- predict(x,confint=confint)
              if (cumulative) pred$deaths <- cumsum(pred$deaths)
              if (show.pred) {
                  with(pred,lines(date,deaths,lwd=lwd.pred,col=col.pred))
                  if (confint) {
                      with(pred,polygon(c(date,rev(date)),
                                        c(lwr,rev(upr)),
                                        border=NA,
                                        col=col.ci))
                  }
              }
              ## display metadata for fit on upper right of plot
              legpos <- if (cumulative || window.only)
                "bottomright" else "topright"
              gr <- growthRate(x)
              if (!add && show.legend) legend(legpos,xpd=NA,bty="n",
                     legend=c(x@model@name,x@loglik@name))
              if (!add && show.gof) legend("topleft",xpd=NA,bty="n",
                     legend=c(sprintf("dates: %s - %s",min(x@date),max(x@date)),
                       sprintf("gof = %g",gof(x)),
                       sprintf("window = %d:%d",x@first,x@last),
                       sprintf("peak = %d",x@peak),
                       sprintf("r = %.3g (%.3g,%.3g)",gr[1],gr[2],gr[3])))
              if (!add && show.totals) {
                ylo <- par("usr")[3]
                yhi <- par("usr")[4]
                yleg <- ylo + 0.4*(yhi-ylo)
                legend(x=par("usr")[1],y=yleg,bty="n",title="Total deaths:",
                       legend=c(sprintf("%d in all",as.integer(sum(x@deaths[seq_along(x@time)],na.rm=TRUE))),
                         sprintf("%d in window",as.integer(sum(x@deaths[x@window])))))
              }
            })

#' predict method
#' @param object \code{\link{epigrowthfit}} object
#' @param tvec time vector
#' @param newpars modified model parameters for prediction
#' @param confint return confidence intervals?
#' @param level confidence interval
#' @param seed random-number seed (for importance sampling CIs)
#' @param nsim number of importance samples for CIs
#' @return a data frame with columns \code{time} and \code{deaths}
#' and (if confidence intervals requested) \code{lwr} and \code{upr}
#' @importFrom bbmle coef pop_pred_samp
#' @importFrom stats quantile
#' @importFrom Hmisc wtd.quantile
#' @export
#' 
setMethod("predict",
          signature(object="epigrowthfit"),
          definition = function(object,
                                tvec=NULL,
                                newpars=NULL,
                                confint=FALSE,
                                seed=101,
                                nsim=500,
                                level=0.95) {
    ## FIXME: redundant with "fitted" method?
    if (!is.null(newpars)) {
        pars <- newpars
    } else {
        pars <- bbmle::coef(object@mle2)
    }
    t_internal <-  object@time[object@window]
    if (is.null(tvec)) {
        tvec <- t_internal
    }
    tvec_off <- tvec-t_internal[1]
    if (confint) {
        if (!is.null(seed)) set.seed(seed)
        simtraj <- matrix(NA,nrow=length(tvec),ncol=nsim)
        simpars <- pop_pred_samp(object@mle2,n=nsim,PDify=TRUE,return_wts=TRUE)
        sim_OK <- !all(is.na(simpars))
        if (sim_OK) {
            for (i in 1:nrow(simtraj)) {
                simtraj[,i] <- predict(object,tvec=tvec,
                                       newpars=simpars[i,],
                                       confint=FALSE)[["deaths"]]
            }
            bad_pars <- is.na(simpars[,"wts"])
            simtraj <- simtraj[,!bad_pars]
            wts <- na.omit(simpars[,"wts"])
            envelope <- try(t(apply(simtraj,1,Hmisc::wtd.quantile,
                                weights=wts,
                                normwt=TRUE,
                                probs=c((1-level)/2,(1+level)/2),
                                na.rm=TRUE)),
                            silent=TRUE)
            if (inherits(envelope,"try-error")) {
                envelope <- matrix(NA,ncol=2,nrow=length(tvec))
            }
        }
    } ## if confint
    res <- data.frame(time=tvec,
                      date=safe_date(tvec),
                      deaths=object@mean(pars,tvec_off))
    if (confint && sim_OK) {
        res <- data.frame(res,lwr=envelope[,1],upr=envelope[,2])
    }
    return(res)
})

##' @export
setMethod("coef",
          "epigrowthfit",
          definition=function(object, scale=c("natural","unconstrained"), ...) {
              scale <- match.arg(scale)
              c0 <- object@mle2@coef
              if (scale=="unconstrained") return(c0)
              c1 <- unlist(transformPar(object@model,c0,inverse=FALSE))
              c2 <- unlist(transformPar(object@loglik,c1,inverse=FALSE))
              ## undo x0 scaling step
              if ("K" %in% names(c2)) {
                  c2[["x0"]] <- c2[["x0"]]*c2[["K"]]
              }
              return(c2)
})

##' @export
setMethod("logLik",
          "epigrowthfit",
          function (object, ...) {
    logLik(object@mle2, ...)
})
