#' Data in epigrowthfitTMB
#'
#' @description
#' Epidemic time series data to which growth rates can be fit
#' using the \pkg{epigrowthfitTMB} package.
#'
#' @details
#' Here is a list of available data sets with links to the
#' relevant documentation:
#'
#' \describe{
#'   \item{[`canadacovid`][canadacovid]}{
#'     Daily confirmations of COVID-19 in Canadian provinces and territories
#'     from February 14, 2020 to June 21, 2020.
#'   }
#'   \item{[`husting`][husting]}{
#'     Weekly counts of wills probated in the Court of Husting
#'     during four plague epidemics in 14th century London.
#'   }
#'   \item{[`canterbury`][canterbury]}{
#'     Weekly counts of wills probated in the Prerogative Court of Canterbury
#'     during 24 plague epidemics in 16th and 17th century London.
#'   }
#'   \item{[`londonparish`][londonparish]}{
#'     Weekly counts of burials listed in extant parish registers
#'     during 24 plague epidemics in 16th and 17th century London.
#'   }
#'   \item{[`londonbills`][londonbills]}{
#'     Weekly counts of plague deaths recorded in the London Bills of Mortality
#'     during 24 plague epidemics in 16th and 17th century London.
#'   }
#' }
#'
#' @keywords internal
#' @name epigrowthfitTMB-data
NULL

#' COVID-19 in Canadian provinces and territories
#'
#' @description
#' Daily confirmations of COVID-19 in Canadian provinces and territories
#' from February 14, 2020 to June 21, 2020.
#'
#' @format
#' A data frame with 1677 rows and 3 variables:
#'
#' \describe{
#'   \item{`date`}{Date formatted as YYYY-MM-DD, daily from
#'     February 14, 2020 to June 21, 2020 within each level of `province`.
#'   }
#'   \item{`time`}{Decimal date, equal to year plus fraction
#'     of year. The fraction is 0 on January 1, 1/365 or 1/366
#'     on January 2, and so on.
#'   }
#'   \item{`province`}{A factor variable with 13 levels
#'     specifying a Canadian province or territory.
#'   }
#'   \item{`new_confirmations`}{Number of new confirmations of COVID-19.
#'     There are 37 missing values (at least one for each level of `province`).
#'   }
#' }
#'
#' @usage data(canadacovid)
#' @examples
#' data(canadacovid)
#' subset(canadacovid, province == "ON")
#'
#' @name canadacovid
"canadacovid"

#' Husting wills
#'
#' @description
#' Weekly counts of wills probated in the Court of Husting
#' during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A data frame with 246 rows and 5 variables:
#'
#' \describe{
#'   \item{`date`}{Date formatted as YYYY-MM-DD, weekly within
#'     each level of `outbreak` with longer gaps between outbreaks.
#'   }
#'   \item{`time`}{Decimal date, equal to year plus fraction
#'     of year. The fraction is 0 on January 1 and either 1/365
#'     or 1/366 on January 2.
#'   }
#'   \item{`wills`}{Count of wills written. Within a given level
#'     of `outbreak`, `wills[i]` is the number of wills written
#'     between `date[i-1]+1` and `date[i]`.
#'   }
#'   \item{`outbreak`}{An ordered factor variable with 4 levels,
#'     partitioning the time series into distinct plague outbreaks.
#'     Levels are named roughly according to the years in which
#'     the outbreaks took place: 1348, 1361, 1368, and 1375.
#'   }
#'   \item{`severity`}{An ordered factor variable with levels `"minor"`
#'     and `"major"`, indicating outbreak severity. All 14th century
#'     outbreaks are classified as major, only because the data are too
#'     sparse to distinguish between minor and major outbreaks
#'     (in contrast with wills probated in the
#'     [Prerogative Court of Canterbury][canterbury]).
#'   }
#' }
#'
#' @usage data(husting)
#' @examples
#' data(husting)
#' subset(husting, outbreak == 1375)
#' @name husting
"husting"

#' Canterbury wills
#'
#' @description
#' Weekly counts of wills probated in the Prerogative Court of Canterbury
#' during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A data frame with 1343 rows and 5 variables:
#'
#' \describe{
#'   \item{`date`}{Date formatted as YYYY-MM-DD, weekly within
#'     each level of `outbreak` with longer gaps between outbreaks.
#'   }
#'   \item{`time`}{Decimal date, equal to year plus fraction
#'     of year. The fraction is 0 on January 1 and either 1/365
#'     or 1/366 on January 2.
#'   }
#'   \item{`wills`}{Count of wills written. Within a given level
#'     of `outbreak`, `wills[i]` is the number of wills written
#'     between `date[i-1]+1` and `date[i]`.
#'   }
#'   \item{`outbreak`}{An ordered factor variable with 24 levels,
#'     partitioning the time series into distinct plague outbreaks.
#'     Levels are named roughly according to the years in which
#'     the outbreaks took place: 1563, 1578, ..., 1647, and 1665.
#'   }
#'   \item{`severity`}{An ordered factor variable with levels
#'     `"minor"` and `"major"`, indicating outbreak severity.
#'     An outbreak is classified as major if and only if plague
#'     deaths per week per 1000 individuals exceeded 5 at least once.
#'   }
#' }
#'
#' @usage data(canterbury)
#' @examples
#' data(canterbury)
#' subset(canterbury, outbreak == 1665) # Great Plague of London
#' @name canterbury
"canterbury"

#' Burials in London parish registers
#'
#' @description
#' Weekly counts of burials listed in extant parish registers
#' during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A data frame with 1320 rows and 5 variables:
#'
#' \describe{
#'   \item{`date`}{Date formatted as YYYY-MM-DD, weekly within
#'     each level of `outbreak` with longer gaps between outbreaks.
#'   }
#'   \item{`time`}{Decimal date, equal to year plus fraction
#'     of year. The fraction is 0 on January 1 and either 1/365
#'     or 1/366 on January 2.
#'   }
#'   \item{`burials`}{Count of burials. Within a given level
#'     of `outbreak`, `burials[i]` is the number of burials
#'     registered between `date[i-1]+1` and `date[i]`.
#'   }
#'   \item{`outbreak`}{An ordered factor variable with 24 levels,
#'     partitioning the time series into distinct plague outbreaks.
#'     Levels are named roughly according to the years in which
#'     the outbreaks took place: 1563, 1578, ..., 1647, and 1665.
#'   }
#'   \item{`severity`}{An ordered factor variable with levels
#'     `"minor"` and `"major"`, indicating outbreak severity.
#'     An outbreak is classified as major if and only if plague
#'     deaths per week per 1000 individuals exceeded 5 at least once.
#'   }
#' }
#'
#' @usage data(londonparish)
#' @examples
#' data(londonparish)
#' subset(londonparish, outbreak == 1665) # Great Plague of London
#' @name londonparish
"londonparish"

#' Plague deaths in the London Bills of Mortality
#'
#' @description
#' Weekly counts of plague deaths recorded in the London Bills of Mortality
#' during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A data frame with 1329 rows and 7 variables:
#'
#' \describe{
#'   \item{`date`}{Date formatted as YYYY-MM-DD, weekly within
#'     each level of `outbreak` with longer gaps between outbreaks.
#'   }
#'   \item{`time`}{Decimal date, equal to year plus fraction
#'     of year. The fraction is 0 on January 1 and either 1/365
#'     or 1/366 on January 2.
#'   }
#'   \item{`all_causes_deaths`}{Count of deaths due to all causes.
#'     Within a given level of `outbreak`, `all_cause_deaths[i]`
#'     is the number of all causes deaths between `date[i-1]+1`
#'     and `date[i]`. There are 30 missing values, all during
#'     `outbreak = 1563`.
#'   }
#'   \item{`plague_deaths`}{Count of deaths due to plague.
#'     Within a given level of `outbreak`, `plague_deaths[i]`
#'     is the number of plague deaths between `date[i-1]+1`
#'     and `date[i]`.
#'   }
#'   \item{`outbreak`}{An ordered factor variable with 24 levels,
#'     partitioning the time series into distinct plague outbreaks.
#'     Levels are named roughly according to the years in which
#'     the outbreaks took place: 1563, 1578, ..., 1647, and 1665.
#'   }
#'   \item{`severity`}{An ordered factor variable with levels
#'     `"minor"` and `"major"`, indicating outbreak severity.
#'     An outbreak is classified as major if and only if plague
#'     deaths per week per 1000 individuals exceeded 5 at least once.
#'   }
#'   \item{`population`}{Estimated London population size.}
#' }
#'
#' @usage data(londonbills)
#' @examples
#' data(londonbills)
#' subset(londonbills, outbreak == 1665) # Great Plague of London
#' @name londonbills
"londonbills"

