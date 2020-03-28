#' Get range
#' @param ep endpoints (in decimal time)
#' @param data data frame
#' @param time_col name of time variable
#' @param maxlen max window length (no longer used)
#' @export
get_range <- function(ep,data,maxlen=Inf,time_col="time") {
    s <- seq(which.min(abs(data[[time_col]]-ep[1])),
             which.min(abs(data[[time_col]]-ep[2])))
    if (length(s)>maxlen) {
        s <- s[1:maxlen]
    }
    return(s)
}

#' Get baseline
#' 
#' add an estimated deaths column that is the total deaths minus a baseline estimated by subtracting the average mortality in surrounding years for exactly the part of the year we are using for our fits
#' @param data data frame including at least columns \code{time} and \code{in_col}
#' @param epid_defs data frame including outbreak.year, start, end
#' @param in_col name of deaths column in input
#' @param out_col name of deaths column in output
#' @param time_col name of numeric time column in input
#' @param per_wid number of years on either side of epidemic to include baseline
#' @param min_val minimum value (floor) for baseline
#' @param include_all include all values, not just epidemic-year data (i.e. time series of counts between epidemics) in the returned data frame?
#' @param do_plot produce diagnostic plots?
#' @param y_min minimum y-value for plotting
#' @return a data frame of similar structure to the one that was input as \code{data}, with the following modifications:
#' \itemize{
#' \item new columns \code{prev<i>} and \code{next<i>} (where <i> runs from 1 to \code{perwid} representing the counts in previous and succeeding years; these are useful for diagnostic plots
#' \item a column labelled \code{out_col} (the difference between \code{in_col}
#' and the computed baseline)
#' \item an attribute \code{baseline}: a named numeric vector with the baseline value for each epidemic in \code{epid_defs} (names correspond to the \code{outbreak.year} variable in \code{epid_defs}
#' }
#' @details for the method to work properly, the data should include values for the periods before and after each outbreak year included in \code{epidemic_defs} (depends on \code{per_wid}: +/- 2 years by default)
#' @export
#' @examples
#' \dontrun{
#' library(dplyr)
#' ctmp <- canterbury_wills_individual %>%
#'    get_date_count("2 weeks") %>%
#'    get_baseline(epidemic_defs,in_col="nwills",
#'                 do_plot=FALSE,include_all=TRUE)
#' }
get_baseline <- function(data,epid_defs,
                         in_col="total",
                         out_col="plague.deaths",
                         time_col="time",
                         per_wid=2,
                         min_val=0,
                         include_all=FALSE,
                         do_plot=FALSE,
                         y_min=0.1) {

    if (!in_col %in% names(data)) {
        stop("input column not found in data")
    }
    majors <- epid_defs$outbreak.year
    
    ## make empty list of appropriate length with appropriate names:
    dat0 <- setNames(vector("list",length(majors)),majors)

    ## generate baseline, make corrections, and plot diagnostics:

    for (y in majors) {
        ## if (y==1665) browser()
        ## should be: 1664.989 1666.962
        endpoints <- epid_defs[epid_defs$outbreak.year==y,c("start","end")]
        endpoints <- unlist(endpoints)
        ## endpoints <- unlist(subset(epid_defs,outbreak.year==y,
        ##            select=c(start,end)))
        yn <- as.character(y)
        tvals <- get_range(endpoints,data, time_col=time_col)
        if (length(tvals)>1) { ## any data for this epidemic
            dat0[[yn]] <- list(data=data.frame(obs=seq_along(tvals),data[tvals,]))
            ## tot_death = total number of deaths in the baseline period
            ## tot_n = total number of observations in the baseline period
            tot_death <- tot_n <- 0
            ## years to include in baseline calculation
            v <- seq(-per_wid,per_wid)
            v <- v[v!=0]

            all_offsets <- character(0)
            for (i in seq_along(v)) {
                if (v[i]<0) {
                    dirstr <- "prev"
                    tvals <- get_range(endpoints[1]+v[i]+c(-1,0),data,
                                       time_col = time_col)
                } else {
                    dirstr <- "next"
                    tvals <- get_range(endpoints[2]+v[i]+c(0,1),data,
                                       time_col = time_col)
                }
                cur_dat <- data.frame(obs=seq_along(tvals),data[tvals,in_col])
                oname <- paste0(dirstr,abs(v[i]))
                names(cur_dat)[2] <- oname
                all_offsets <- c(all_offsets,oname)
                ## need to account for NAs in both number of obs & total deaths
                tot_n <- tot_n + sum(!is.na(cur_dat[[2]]))
                tot_death <- tot_death+sum(cur_dat[[2]],na.rm=TRUE)
                ## FIXME: merge might not match exactly
                ##  (but since we're only using the prev/next values for an
                ##   overall baseline, not an exact match, probably doesn't
                ##   matter that much)
                dat0[[yn]][["data"]] <-
                         merge(dat0[[yn]][["data"]],cur_dat,all.x=TRUE)
            }
            dat0[[yn]]$baseline0 <- tot_death/tot_n
            dat0[[yn]]$baseline <- 
                bd <- round(tot_death/tot_n) ## numbers must be integers for Poisson or NegBin
            if (do_plot) {
                m <- dat0[[yn]]$data[,c(in_col,all_offsets)]
                y_min <- min(min(m[m>0]),bd[bd>0],y_min)
                mlab <- y
                if (bd==0) mlab <- paste(mlab,"(baseline==0)")
                ## suppress log-zero warnings
                suppressWarnings(matplot(m,
                                         xlab="time",
                                         main=mlab,
                                         log="y",
                                         type="b",
                                         col=1:(1+length(v)),
                                         lty=1,pch=1,
                                         ylim=c(y_min,max(m))))
                abline(h=bd,lwd=2,lty=2)
            } ## if do_plot
        } ## if any data
    } ## loop over outbreak years
    ## browser()
    
    ## strip missing values
    dat0 <- dat0[!sapply(dat0,is.null)]

    for (yn in names(dat0)) {
        dat0[[yn]]$data[[out_col]] <-
            with(dat0[[yn]],pmax(min_val,data[[in_col]]-baseline))
        dat0[[yn]]$data$outbreak.year <- as.numeric(yn)
    }
    r <- do.call(rbind,lapply(dat0,"[[","data"))
    r <- r[names(r) != "obs"] ## exclude temporary obs variable
    attr(r,"baseline") <- sapply(dat0,"[[","baseline0")

    if (include_all) {

        ## retrieve data that's *not* in the data file
        other_data <- data[! (data$date %in% r$date), ]

        ## dplyr::bind_rows would be easier, but
        ##  want to avoid unnecessary dependency
        other_data[c("prev2","prev1",
                     "next1","next2",out_col,"outbreak.year")] <- NA

        r <- rbind(r,other_data)
        r <- r[order(r$date),]
    }
    
    rownames(r) <- NULL
    return(r)
}
