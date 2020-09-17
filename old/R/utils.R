##' compute R0 given r and generation interval distribution
##' @param obj instantaneous exponential growth rate
##' @param generation.interval the generation interval distribution (days)
##' @docType methods
##' @export
setMethod(
  "R0", 
  "numeric",
  definition = function(obj, generation.interval=generation.interval.plague()) {
    if (length(obj) == 1) { 
      t = seq(0, by=1, along.with=generation.interval)/365
      1 / sum(exp(-obj*t)*generation.interval)
    } else sapply(obj, R0)
  }
)

##' @export
setMethod(
  "R0", 
  "matrix",
  definition = function(obj, ...) {
      obj[] <- R0(c(obj))
      return(obj)
  }
)

#' Generation interval distribution for pneumonic plague
#'
#' @description
#' compute generation interval for plague by convolving
#' (built-in) latent period distribution with reversed
#' infectious period distribution
#'
#' @references
#' \insertRef{Sven07}{epigrowthfit}
#' @importFrom stats convolve
#' @export
generation.interval.plague <- function() {
  convolve(epigrowthfit::latent.period$frequency,
           rev(epigrowthfit::infectious.period$frequency),
           type="open")
}

#' compute doubling time from exponential growth rate
#'
#' @description
#' The doubling time in days is \code{log(2)/obj*365},
#' where \code{obj} is the exponential growth rate in years.
#' @param obj growth rate in units 1/year
#' @export
setMethod(
  "doublingTime",
  "numeric",
  definition = function(obj) { log(2)/obj*365 }
)

#' expected final size of an epidemic
#'
#' @description
#' For a given R0, return the proportion of the population
#' expected to be infected if an epidemic is seeded and
#' all individuals are initially susceptible.  The standard
#' final size formula (Kermack and McKendrick 1927,
#' Ma and Earn 2006) is used.
#' @references
#' \insertRef{KermMcKe27}{epigrowthfit}
#'
#' \insertRef{MaEarn06}{epigrowthfit}
#' 
#' @param R0 basic reproduction number
#' @export
#' @importFrom emdbook lambertW
## FIXME: why does finalsize(1.13796303983581) warn?
## (finalsize(c(1.135,1.140)) is fine ...
finalsize <- function(R0) {
    ## don't know 
    fs <- 1 + 1/R0 * suppressWarnings(lambertW(-R0*exp(-R0), maxiter=1000))
    return(fs)
}

#########################################
## plotPlague and associated functions ##
#########################################

## used in plotPlague
plotPlague.xaxis.labels <- function(data) {
  is.monthly <- (data$aggregation[1]=="monthly")
  bars.per.year <- if (is.monthly) 12 else 52
  mnames <- month.abb
  ## if we are plotting more than two years then only label
  ## every second month:
  if (length(unique(data$year)) >= 3) mnames[which((1:12)%%2 == 0)] <- ""
  if (is.monthly) {
    ## every month
    labels <- mnames[data$month]
  } else {
    ## first week of each month
    labels <- with(data,ifelse(day<8,mnames[month],""))
  }
  return(labels)
}

## used in fancy.axis.Date
add.yearlabs.to.axis.Date <- function(df, col.year="black") {
    with(df,{
        nyears <- diff(range(year)) + 1
        if (nyears > 1 && nyears < 6) { # indicate start of each year
            yearlabs <- as.character(unique(year))
            mtext(yearlabs,side=1,line=2,cex=1,col=col.year,
                  at=as.Date(sprintf("%d-01-01",unique(year))))
        }
    })
}

## used in plotPlague
add.yearlabs.to.barplot <- function(data,epidemic,bar.midpts) {
  nbars <- nrow(data)
  with(data,{
    bars.per.year <- if (aggregation[1]=="monthly") 12 else 52
    if (nbars > bars.per.year) { # indicate start of each year
      ## indices of January or first week of January:
      iJan <- if (aggregation[1]=="monthly")
        which(data$month == 1) else which(data$month == 1 & data$day <= 7)
      if (length(epidemic) == 1) {
        ## this is a single epidemic spanning multiple years,
        ## so mark the start of each year with the year:
        yearlabs <- as.character(epidemic + 0:(length(iJan)-1))
      } else {
        ## multiple epidemics are plotted, so just indicate
        ## Januaries with a vertical bar:
        yearlabs <- "|"
      }
      mtext(yearlabs,side=1,line=3,at=bar.midpts[iJan],cex=1)
    }
  })
}

## used in plotPlague
## http://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string
## N.B. this is NOT vectorized
## as.character is essential because the arguments to which this is applied are factors
simpleCap <- function(x) {
  ##s <- strsplit(x, " ")[[1]]
  s <- strsplit(as.character(x), " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}
## Ben's clever solution that doesn't involve understanding the function:
simpleCap <- Vectorize(simpleCap)

## used in plotPlague
## opposite of %in%
## http://stackoverflow.com/questions/5831794/opposite-of-in
'%!in%' <- function(x,y)!('%in%'(x,y))

## used in plotPlague
plotPlague.toprightlabel <- function(epidemic,toprightlabel) {
  if (is.null(epidemic)) {
    toprightlabel <- "ALL"
  } else if (is.null(toprightlabel)) {
    if (length(epidemic)==1) {
      toprightlabel <- as.character(epidemic)
    } else {
      toprightlabel <- sprintf("%d-%d",min(epidemic),max(epidemic))
    }
  }
  return(toprightlabel)
}

#' convert dates to a data frame of year, month, and day
#' @param datetxt dates as character
#' @export date_to_ymd
date_to_ymd <- function(datetxt) {
  df <- data.frame(year = as.numeric(format(datetxt, format = "%Y")),
                   month = as.numeric(format(datetxt, format = "%m")),
                   day = as.numeric(format(datetxt, format = "%d")))
  return(df)
}

#' guess aggregation level of incidence data
#' @param data data frame containing a \code{day} column
#' @return a text label
#' @export
aggregation.level <- function(data) {
  if (!is.null(aggattr <- attr(data,"aggregation"))) return(aggattr)
  ## FIXME: could probably be
  ##   dt <- diff(data$day)[1]
  ##   lab <- switch(dt,"1"=daily,"7"="weekly",...,
  ##             sprintf("%d-daily",dt))
  ##   return(lab)
  with(data,{
    if (diff(day)[1] == 21) return("tri-weekly")
    if (diff(day)[1] == 14) return("bi-weekly")
    if (diff(day)[1] == 7) return("weekly")
    if (diff(day)[1] == 1) return("daily")
    return(sprintf("%d-daily",diff(day)[1]))
  })
}

#' add columns to data frame so it looks like original plague data frame
#' @param data data frame
#' @export
add.plague.data.frame.columns <- function(data) {
  cn <- colnames(data)
  if ("time" %in% cn) {
    if (!("outbreak.year" %in% cn)) data$outbreak.year <- round(data$time)
    if (!("date" %in% cn)) data$date <-
      ##with(data,as.Date(sprintf("%d-%2d-%2d",year,month,day)))
      as.Date(lubridate::date_decimal(data$time))
    if (!("year" %in% cn)) {
      ymdcols <- date_to_ymd(data$date)
      data[names(ymdcols)] <- ymdcols ## avoid losing attributes
    }
    if (!("type" %in% cn)) data$type <- attr(data,"type")
    if (!("place" %in% cn)) data$place <- attr(data,"place")
    if (!("severity" %in% cn)) data$severity <- "?severity?"
    if (!("aggregation" %in% cn)) data$aggregation <- aggregation.level(data)
  }
  if (!("all.cause.deaths" %in% cn)) data$all.cause.deaths <- NA
  return(data)
}

## used in plotPlague
plotPlague.setup <- function(epidemic,data,name.data="plague.deaths") {
  data <- add.plague.data.frame.columns(data)

  plague <- data # FIX: use 'data' rather than 'plague' as the data frame name
  all.years <- unique(plague$outbreak.year)
  if (is.null(epidemic)) {
    epidemic <- all.years
  } else {
    epidemic <- all.years[ which(min(epidemic) <= all.years & all.years <= max(epidemic)) ]
  }
  if (length(epidemic) > 1)
    warning(sprintf("plotPlague.setup: plotting %d epidemics",
                    length(unique(epidemic[which((epidemic %in% plague$outbreak.year))]))))
  if (!all(epidemic %in% plague$outbreak.year)) {
    ## epidemics not found in plague data frame
    ## (DE: I think it is no longer possible to get here)
    missing.epidemics <- paste(epidemic[which(epidemic %!in% plague$outbreak.year)],collapse=",")
    warning("not found in plague data frame: ",
             missing.epidemics)
    epidemic <- unique(epidemic[which(epidemic %in% plague$outbreak.year)])
    warning("plotting the following epidemics: ", paste(epidemic,collapse=","))
  }
  data2 <- plague[plague$outbreak.year %in% epidemic,]
  if (length(unique(data2$aggregation)) > 1)
    warning(sprintf("plotPlague.setup: %d aggregation levels shown together",length(unique(data2$aggregation))))
  if (length(unique(data2$type)) > 1)
    warning(sprintf("plotPlague.setup: %d data types shown together",length(unique(data2$type))))
  return(data2)
}

## used in plotPlague and plot method for epigrowthfit
## df must have year and date columns
#' @export
fancy.axis.Date <- function(df, col.year="black") {
    with(df,{
        nyears <- diff(range(year)) + 1
        desired.n <- 12
        min.n <- 6
        axis.line <- NA # the default value of 'line' for axis()
        ## see ?pretty.Date and ?strptime for format:
        ## Note: after doing all this I discovered lubridate::pretty_dates
        ##       which might be better...
        if (nyears > 5 ) {
          nice.format <- "%Y"
          desired.n <- nyears
          min.n <- 5
        } else if (diff(range(time)) > 0.6) {
          nice.format <- "%b"
        } else {
          nice.format <- "%e %b"
        }
        xtickpos <- axis.Date(side=1,
                              at=pretty(date, n=desired.n, min.n=min.n),
                              ##at=lubridate::pretty_dates(date,desired.n),
                              format=nice.format, line=axis.line)
        add.yearlabs.to.axis.Date(df, col.year=col.year)
        ## minor ticks on x axis (complicated because we want them at months):
        ## extra (minor) ticks at each month for multi-year plots:
        if (nyears > 1 && nyears < 8) {
            at <- c()
            for (yy in unique(year)) {
                at <- c(at,as.Date(sprintf("%d-%2d-01",rep(yy,11),2:12)))
            }
            ## rug() gives clipping warnings
            suppressWarnings(rug(x=at, ticksize=-0.02, xpd=FALSE))
            at.as.Date <- as.Date(at, origin=as.Date("1970-01-01"))
            return(invisible(at.as.Date))
        } else {
          warning("fancy.axis.Date: returning xtickpos rather than at")
          return(invisible(xtickpos))
        }
    })
}

#' plot plague epidemics
#'
#' @description
#' A plot of one or more plague epidemics available
#' in the \code{\link{plague}} data frame.
#'
#' @details
#' If multiple epidemic years are given then all will be plotted.
#'
#' WARNINGS:
#' \itemize{
#' \item With \code{plot.type="barplot"}, temporal gaps between
#'       distinct epidemics will not appear.
#' \item If multiple epidemics are plotted and they include
#'       different levels of aggregation (e.g., weekly, monthly)
#'       the wrong impression will be given for the magnitude
#'       of more highly aggregated epidemics.
#' }
#' 
#' @param epidemic epidemic year(s); all available years by default
#' @param data data frame containing epidemic time series
#' @param name.data the column name of the variable to be plotted (\code{"plague.deaths"} by default)
#' @param toprightlabel epidemic name by default
#' @param cex.toprightlabel character expansion factor for \code{toprightlabel}
#' @param col.toprightlabel font colour for \code{toprightlabel}
#' @param col.data fill colour for points in \code{\link{plot}} or
#'        bars in \code{\link{barplot}}
#' @param cex.data character expansion factor for plague data points
#' @param show.place logical: if \code{TRUE} then indicate place in y axis label
#' @param cumulative logical: if \code{TRUE} then show cumulative rather than raw data
#' @param show.acm logical: if \code{TRUE} then show all cause mortality
#' @param show.legend logical: if \code{TRUE} then show a legend indicating which curve is which (relevant only if \code{show.acm==TRUE})
#' @param col.acm colour of all cause mortality data if shown
#' @param cex.acm character expansion factor for all cause mortality points
#' @param plot.type character string: either "barplot" or "timeplot"
#' @param col.year colour of year labels on x axis (for multi-year plots)
#' @param log log scale ("", "x", "y", or "xy", as in \code{\link{plot.default}})
#' @param add (logical) add to an existing plot?
#' @param ... further parameters passed to \code{\link{plot}}, \code{\link{barplot}}, or \code{\link{points}}
#' @return A \code{\link{data.frame}} that is the subset of the \code{\link{plague}} data
#' frame used to make the plot, with an additional \code{date} column of class \code{\link{Date}}.
#' @examples
#' ## plot the Black Death in London
#' plotPlague(1348)
#' plotPlague(1348,plot.type="barplot")
#' ## plot a sequence minor plague epidemics in London
#' plotPlague(1578:1582)
#' plotPlague(1578:1582,plot.type="barplot")
#' ## plot the Great Plague of London
#' plotPlague(1665)
#' @importFrom graphics points curve legend abline
#' @importFrom Hmisc minor.tick
#' @export
plotPlague <- function(epidemic=NULL, # FIX: should probably call this outbreak.year as in plague data frame
                       data=epigrowthfit::plague,
                       name.data="plague.deaths",
                       toprightlabel=NULL,
                       cex.toprightlabel=2,
                       col.toprightlabel="black",
                       col.data="red",
                       cex.data=1,
                       show.place=TRUE,
                       cumulative=FALSE,
                       show.acm=FALSE, # all cause mortality
                       show.legend=show.acm,
                       col.acm="grey",
                       cex.acm=0.5,
                       plot.type="timeplot", # or "barplot"
                       col.year="blue",
                       add=FALSE,
                       log="",
                       ...) {
  toprightlabel <- plotPlague.toprightlabel(epidemic,toprightlabel)
  data <- plotPlague.setup(epidemic,data,name.data)
  vals.data <- data[,name.data]
  with(data,{
    yvals <- if (cumulative) cumsum(vals.data) else vals.data
    acm <- if (cumulative) cumsum(all.cause.deaths) else all.cause.deaths
    if (plot.type=="barplot") {
      if (show.acm) warning("plotPlague: no barplot option for acm")
      if (add) warning("plotPlague: cannot add to existing barplot")
      xaxis.labs <- if (any(toprightlabel=="ALL")) NA else plotPlague.xaxis.labels(data)
      bar.midpts <- barplot(yvals,col=col.data,
                            names.arg=xaxis.labs,
                            las=2,
                            log=log,...)
    } else {
      if (!add) {
        ymax <- if (show.acm) max(acm,yvals,na.rm=TRUE) else max(yvals,na.rm=TRUE)
        ymin <- if (grepl("y",log)) 1 else 0
        plot(range(date),c(ymin,ymax),type="n",
             las=1,xlab="",ylab="",
             xaxt="n", # xaxis is made using fancy.axis.Date below
             log=log,...)
        ## grey dotted lines at start of each year:
        abline(v=unique(year),col="grey",lty="dotted")
        ## minor ticks on y axis:
        Hmisc::minor.tick(nx=1,ny=2) # one minor tick between adjacent major ticks
        ## x axis:
        fancy.axis.Date(data, col.year=col.year)
      }
      if (show.acm) points(date,acm,pch=21,bg=col.acm,cex=cex.acm,type="o")
      points(date,yvals,pch=21,bg=col.data,cex=cex.data,type="o",...)
    }
    if (!add) {
      legend("topright",
             legend=toprightlabel,cex=cex.toprightlabel,bty="n",
             text.col="black",xpd=NA)
      ylab <- paste(simpleCap(aggregation[1]), simpleCap(type[1]))
      if (show.place) ylab <- paste(ylab, "in", place[1])
      if (cumulative) ylab <- paste("Cumulative", ylab)
      if (plot.type!="barplot") title(ylab=ylab)
      mtext(ylab,side=3,line=0.5,at=0,adj=0.25,cex=1.5,xpd=NA)
      if (plot.type=="barplot") add.yearlabs.to.barplot(data,epidemic,bar.midpts)
      if (show.legend) {
        if (!show.acm) warning("plotPlague: showing legend without acm???")
        legend("topleft",legend=c("all causes","plague"),bty="n",
               pt.bg=c(col.acm,col.data),pt.cex=c(cex.acm,cex.data),
               lty=1,pch=21)
      }
    }
  })
  return(invisible(data))
}

#' tabulate numbers of wills & testaments
#' @param dd data frame
#' @param unit time unit for aggregation
#' @param date_col name of date column
#' @param count_col name of aggregated column
#' @param start_weekday starting weekday for aggregation
#' @param start_date starting date for aggregation
#' @inheritParams aggregate_wills
#' @importFrom lubridate decimal_date
#' @export
get_date_count <- function(data,
                           unit="1 days",
                           date_col="date",
                           count_col="nwills",
                           start_weekday=NA,
                           start_date=NA) {
    ##
    date <- data[[date_col]]
    if (!is.na(start_weekday) || !is.na(start_date)) {
        if (!is.na(start_weekday)) {
            first_day <- date[1]
            ## compute weekday (Monday==1. Sunday==7)
            first_weekday <- as.numeric(format(first_day,"%u"))
            ## figure out how far back we have to go to start bins on a specified weekday
            date_offset <- (start_weekday-first_weekday) %% 7
            start_date <- first_day-date_offset
        }
        date_breaks <- seq(start_date, tail(date,1),
                           by = unit)
        data$fdate <- cut(date, breaks=date_breaks)
    } else {
        data$fdate <- cut(date, breaks=unit)
    }
    data2 <- setNames(as.data.frame(table(data$fdate)),c(date_col,count_col))
    data2[[date_col]] <- as.Date(as.character(data2[[date_col]]))
    data2$time <- lubridate::decimal_date(data2[[date_col]])
    return(data2)
}

## https://stackoverflow.com/questions/23724815/roxygen2-issue-with-exporting-print-method
#' @method summary epidemic.data
#' @export
summary.epidemic.data <- function(object, ...) {
    sub.names <- c("nwills","total","plague.deaths")
    sub.names <- intersect(sub.names,names(object))
    sub.object <- object[,sub.names,drop=FALSE]
    ss <- summary.data.frame(sub.object, ...)
    for (a in c("type","place","aggregation","source")) {
        attr(ss,a) <- attr(object,a)
    }
    attr(ss,"timerange") <- range(object$time,na.rm=TRUE)
    if ("date" %in% names(object)) {
        attr(ss,"daterange") <- range(object$date,na.rm=TRUE)
    }
    attr(ss,"outbreak.year") <- unique(object$outbreak.year)
    attr(ss,"nrow") <- nrow(object)
    class(ss) <- c("summary.epidemic.data","table")
    return(ss)
}

#' @method print summary.epidemic.data
#' @export
print.summary.epidemic.data <- function(x,...) {
    cat("data type:   ",attr(x,"type"),"\n")
    cat("place:       ",attr(x,"place"),"\n")
    cat("aggregation: ",attr(x,"aggregation"),"\n")
    cat("source:      ",attr(x,"source"),"\n")
    cat("time range:  ", attr(x,"timerange"), "\n")
    if (!is.null(dr <- attr(x,"daterange")))
        cat("date range:  ", sprintf("%s", dr), "\n")
    if (!is.null(oy <- attr(x,"outbreak.year")))
        cat("outbreak years: ", oy, "\n")
    cat("number of records: ", attr(x,"nrow"), "\n")
    ## for individual records, there will be no nwills or plague.deaths column:
    if (ncol(x) > 0) {
        NextMethod(x) # print the usual summary of a data frame
    }
}

#' @export
subset.epidemic.data <- function(x, ...) {
    aa <- attributes(x)
    res <- NextMethod(x,...)
    for (a in c("type","place","aggregation","source")) {
        attr(res,a) <- aa[[a]]
    }
    attr(res,"baseline") <-
        aa$baseline[as.character(unique(res$outbreak.year))]
    class(res) <- c("epidemic.data","data.frame")
    return(res)
}
#' aggregate individual wills into a time series
#'
#' @description
#' Aggregate line-listed wills into a time series
#' of a specified temporal resolution.
#'
#' @param wills_individual line-listed file containing individual wills
#' @param aggregation aggregation level
#' @param aggname text to display when referring to the aggregation level
#' @param start_weekday numeric value of weekday (1-7, Monday==1, but note that weekday assignment is wrong (but consistent) pre-Gregorian calendar anyway) on which to start binning (default NA means start on first day of dataset, whatever weekday that is)
#' @return A \code{\link{data.frame}} that is also of class \code{epidemic.data}
#' @examples
#' husting_wills_3weeks <- aggregate_wills(husting_wills_individual,"3 weeks","tri-weekly")
#' plotPlague(data=husting_wills_3weeks)
#'
#' @export
aggregate_wills <- function(wills_individual,
                            aggregation="1 weeks",
                            aggname=aggregation,
                            include_all=FALSE,
                            start_weekday=NA) {

    w1 <- get_date_count(wills_individual,
                         aggregation,start_weekday = start_weekday)
    wills <- get_baseline(w1, epigrowthfit::epidemic_defs,
                          in_col="nwills",do_plot=FALSE,
                          include_all=include_all)
  
    attr(wills,"type") <- "wills"
    attr(wills,"place") <- "London"
    attr(wills,"aggregation") <- aggname
    attr(wills,"source") <- attr(wills_individual, "source")
  
    class(wills) <- c("epidemic.data","data.frame")
    return(wills)
}

#' @export
extract_baseline <- function(data, yr=NULL, unit=c("year","step"),
                             time_col="time") {
    unit <- match.arg(unit)
    ## cases per reporting step ...
    bb <- attr(data,"baseline")
    if (!is.null(yr)) {
        bb <- bb[as.character(yr)]
    }
    if (unit=="year") {
        ## time step ...
        dt <- diff(data[[time_col]])[1]
        ## -> cases/year
        bb <- bb/dt
    }
    return(bb)
}

#' Take data that are already in count form and change time unit
#' Should be faster than the other alternative
#' @export
#' @importFrom zoo zoo as.yearmon
#' @importFrom stats aggregate
aggsum_monthly <- function(year,month,day,deaths) {
    datevec <- as.Date(paste(year,month,day,sep="-"))
    pdz <- zoo(deaths,order.by=datevec)
    return(aggregate(pdz,as.yearmon,sum))
}

#' aggregate count data into a different time unit
#' (don't try increasing the frequency this way!)
#' @param year year
#' @param month month
#' @param day day
#' @param deaths number of counts (of whatever)
#' @param period time period for aggregation
#' @param time decimal date (alternative to year/month/day)
#' @param percorrect rescale to original units
#' @examples
#' agg1 <- with(london_bills,aggsum_monthly(year,month,day,plague.deaths))
#' agg2 <- with(london_bills,aggsum(year,month,day,plague.deaths))
#' agg3 <- with(london_parish_all,aggsum(time=time,deaths=total))
#' plot(deaths~time,data=agg3,type="l",log="y")
#' agg4 <- with(london_parish_all,aggsum(time=time,deaths=total,
#'                                       period="10 days"))
#' plot(deaths~time,data=agg4,type="l",log="y")
#' @export
aggsum <- function(year,month,day,deaths,period="1 months",time=NULL,
                   percorrect=TRUE) {
    if (!is.null(time)) {
        datevec <- safe_date(time)
    } else {
        datevec <- as.Date(paste(year,month,day,sep="-"))
    }
    aggdate <- seq.Date(min(datevec),max(datevec)+30,by=period)
    aggperiod <- cut.Date(datevec,aggdate)
    deaths <- tapply(deaths,list(aggperiod),sum,na.rm=TRUE)
    if (percorrect) {
        aggtab <- table(aggperiod)

        deaths <- deaths*max(aggtab)/aggtab ## ... dangerous!
        ## but not sure what to do:
        ## ceiling(mean(aggtab)) gives weird results
    }
    res <- data.frame(date=aggdate,
                      time=lubridate::decimal_date(aggdate),
                      deaths=c(deaths,NA))
    return(res)
}

##' summarize nested lists of fits
##' 
##' level=0 returns the summary from a single fit
##' level=1 returns the summary from a list of fits
##' level=2 returns the summary from a nested list of fits ...
##' @param object a (possibly nested) list of epidemic fits
##' @param level number of levels to collapse
##' @param keys names of levels
##' @param type "sum": summary; "pred"; predicted values; "deaths"; original data used to fit
##' @param gof_agg aggregation unit for goodness of fit calculations
##' @param verbose print messages?
##' @param confint_pred ??
##' @param add_NBdisp_coef add NA value for NB disp if necessary?
##' @param ... extra arguments (ignored: for generic consistency)
##' @importFrom methods as
##' @method summary fitList
##' @export summary.fitList
summary.fitList <- function(object,level=1,keys=c("outbreak.year","source"),
                            type="sum",gof_agg="7 days", verbose=FALSE,
                            confint_pred=TRUE, add_NBdisp_coef=TRUE, ...) {
    maxLevel <- length(keys)+1
    b <- if (level==0) object else object[[rep(1,level)]]
    template <- as(summary(b),"data.frame")
    template[] <- NA ## replace *contents* of template with NA
    funList <- list()
    sumfun <- function(x) {
        r <- as(summary(x,gof_agg=gof_agg,...),"data.frame")
        ## FIXME/hack: would be cleaner to change bindrows() below
        ## to use dplyr::bindrows or plyr::rbind.fill, but trying
        ## not to add dependencies ...
        if (add_NBdisp_coef && r$distribname[1]=="poisson") {
            r$coefficients.ll.k <- NA
        }
        return(r)
    }
    funList[[1]] <- function(x,type="sum") {
        if (length(x)==0 || inherits(x,"try-error")) {
            if (type=="sum") ss <- template else ss <- NULL
        } else {
            ss <- switch(type,
                         sum=sumfun(x),
                         pred=predict(x,confint=confint_pred),
                         deaths=data.frame(time=x@time,deaths=x@deaths))
        }
        return(ss)
    }
    ## replacement for dplyr::bind_rows
    ## (doesn't handle mismatched columns though ...)
    bindrows <- function(L,.id) {
        nn <- function(x) {
            if (is.null(n <- nrow(x))) 0 else n
        }
        s <- vapply(L,nn,numeric(1))
        if (length(s)==0 || all(s==0)) return(NULL)
        if (any(s==0)) {
            L <- L[s>0]
            s <- s[s>0]
        }
        dd <- data.frame(.xxx=rep(names(L),s),do.call(rbind,L),
                         stringsAsFactors=FALSE)
        names(dd)[1] <- .id
        return(dd)
            
    }
    for (i in 2:maxLevel) {
        funList[[i]] <- function(x,type="sum") {
            ## dplyr::bind_rows(lapply(x,funList[[i-1]],type=type),.id=keys[i-1])                        
            L <- lapply(x,funList[[i-1]],type=type)
            return(bindrows(L,keys[i-1]))
        }
        ## force evaluation of environment
        environment(funList[[i]]) <- list2env(list(i=i,keys=keys))
    }
    ## funList[[1]](fitlist[[1]][[1]],type="pred")
    ## funList[[2]](fitlist[[1]],type="pred")
    ## funList[[3]](fitlist,type="pred")
    return(funList[[level+1]](object,type=type))
}
