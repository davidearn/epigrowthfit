#' Data in the \pkg{epigrowthfit} package
#'
#' @description
#' Epidemic time series to which growth rates can be fit,
#' and a few other data sets.
#'
#' @details
#' Here is a list of available data sets with links to their documentation.
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
#'   \item{[`plague_latent_period`][plague_latent_period]}{
#'     Distributions (observed and fitted) of the latent period
#'     of pneumonic plague.
#'   }
#'   \item{[`plague_infectious_period`][plague_infectious_period]}{
#'     Distributions (observed and fitted) of the infectious period
#'     of pneumonic plague.
#'   }
#' }
#'
#' @keywords internal
#' @name epigrowthfit-data
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
#'   \item{`date`}{\[Date\] Date, daily from February 14, 2020
#'     to June 21, 2020 within each level of `province`.
#'   }
#'   \item{`time`}{\[numeric\] Decimal date, equal to year plus
#'     fraction of year. The fraction is 0 on January 1, 1/365
#'     or 1/366 on January 2, and so on.
#'   }
#'   \item{`province`}{\[factor, 13 levels\] A Canadian province
#'     or territory (postal abbreviation).
#'   }
#'   \item{`new_confirmations`}{\[numeric\] Number of new confirmations
#'     of COVID-19. There are 37 missing values, at least one for each
#'     level of `province`.
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
#'   \item{`date`}{\[Date\] Date, weekly within each level of
#'     `outbreak` with longer gaps between outbreaks.
#'   }
#'   \item{`time`}{\[numeric\] Decimal date, equal to year plus
#'     fraction of year. The fraction is 0 on January 1, 1/365
#'     or 1/366 on January 2, and so on.
#'   }
#'   \item{`wills`}{\[numeric\] Count of wills written. Within
#'     a given level of `outbreak`, `wills[i]` is the number
#'     of wills written between `date[i-1]+1` and `date[i]`.
#'   }
#'   \item{`outbreak`}{\[ordered factor, 4 levels\] Outbreak label,
#'     partitioning the time series into distinct plague outbreaks.
#'     Levels are named roughly according to the years in which
#'     the outbreaks took place: 1348, 1361, 1368, and 1375.
#'   }
#'   \item{`severity`}{\[ordered factor, 2 levels\] Severity label,
#'     classifying outbreaks as `"minor"` or `"major"`. All 14th
#'     century outbreaks are classified as major, only because
#'     the data are too sparse to distinguish between minor and
#'     major outbreaks (in contrast with wills probated in the
#'     [Prerogative Court of Canterbury][canterbury]).
#'   }
#' }
#'
#' @source
#' These data were transcribed from \insertCite{Shar89;textual}{epigrowthfit}.
#'
#' @references
#' \insertRef{Shar89}{epigrowthfit}
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
#'   \item{`date`}{\[Date\] Date, weekly within each level of
#'     `outbreak` with longer gaps between outbreaks.
#'   }
#'   \item{`time`}{\[numeric\] Decimal date, equal to year plus
#'     fraction of year. The fraction is 0 on January 1, 1/365
#'     or 1/366 on January 2, and so on.
#'   }
#'   \item{`wills`}{\[numeric\] Count of wills written. Within
#'     a given level of `outbreak`, `wills[i]` is the number
#'     of wills written between `date[i-1]+1` and `date[i]`.
#'   }
#'   \item{`outbreak`}{\[ordered factor, 24 levels\] Outbreak label,
#'     partitioning the time series into distinct plague outbreaks.
#'     Levels are named roughly according to the years in which
#'     the outbreaks took place: 1563, 1578, ..., 1647, and 1665
#'     (the Great Plague of London).
#'   }
#'   \item{`severity`}{\[ordered factor, 2 levels\] Severity label,
#'     classifying outbreaks as `"minor"` or `"major"`. An outbreak
#'     is classified as major if and only if plague deaths per week
#'     per 1000 individuals exceeded 5 at least once.
#'   }
#' }
#'
#' @source
#' These data were obtained from the National Archives (UK) in 2018
#' using this
#' [research guide](https://www.nationalarchives.gov.uk/help-with-your-research/research-guides/wills-1384-1858/).
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
#'   \item{`date`}{\[Date\] Date, weekly within each level of
#'     `outbreak` with longer gaps between outbreaks.
#'   }
#'   \item{`time`}{\[numeric\] Decimal date, equal to year plus
#'     fraction of year. The fraction is 0 on January 1, 1/365
#'     or 1/366 on January 2, and so on.
#'   }
#'   \item{`burials`}{\[numeric\] Count of burials. Within a given
#'     level of `outbreak`, `burials[i]` is the number of burials
#'     registered between `date[i-1]+1` and `date[i]`.
#'   }
#'   \item{`outbreak`}{\[ordered factor, 24 levels\] Outbreak label,
#'     partitioning the time series into distinct plague outbreaks.
#'     Levels are named roughly according to the years in which
#'     the outbreaks took place: 1563, 1578, ..., 1647, and 1665
#'     (the Great Plague of London).
#'   }
#'   \item{`severity`}{\[ordered factor, 2 levels\] Severity label,
#'     classifying outbreaks as `"minor"` or `"major"`. An outbreak
#'     is classified as major if and only if plague deaths per week
#'     per 1000 individuals exceeded 5 at least once.
#'   }
#' }
#'
#' @source
#' These data were obtained by \insertCite{Cumm+16;textual}{epigrowthfit}
#' with permission from [Ancestry.com](https://www.ancestry.com/).
#' They have been made publicly available [here](http://www.neilcummins.com).
#'
#' @references
#' \insertRef{Cumm+16}{epigrowthfit}
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
#'   \item{`date`}{\[Date\] Date, weekly within each level of
#'     `outbreak` with longer gaps between outbreaks.
#'   }
#'   \item{`time`}{\[numeric\] Decimal date, equal to year plus
#'     fraction of year. The fraction is 0 on January 1, 1/365
#'     or 1/366 on January 2, and so on.
#'   }
#'   \item{`all_causes_deaths`}{\[numeric\] Count of deaths due
#'     to all causes. Within a given level of `outbreak`,
#'     `all_cause_deaths[i]` is the number of all causes deaths
#'     between `date[i-1]+1` and `date[i]`. There are 30 missing
#'     values, all during `outbreak = 1563`.
#'   }
#'   \item{`plague_deaths`}{\[numeric\] Count of deaths due to plague.
#'     Within a given level of `outbreak`, `plague_deaths[i]` is the
#'     number of plague deaths between `date[i-1]+1` and `date[i]`.
#'   }
#'   \item{`outbreak`}{\[ordered factor, 24 levels\] Outbreak label,
#'     partitioning the time series into distinct plague outbreaks.
#'     Levels are named roughly according to the years in which
#'     the outbreaks took place: 1563, 1578, ..., 1647, and 1665
#'     (the Great Plague of London).
#'   }
#'   \item{`severity`}{\[ordered factor, 2 levels\] Severity label,
#'     classifying outbreaks as `"minor"` or `"major"`. An outbreak
#'     is classified as major if and only if plague deaths per week
#'     per 1000 individuals exceeded 5 at least once.
#'   }
#'   \item{`population`}{\[numeric\] Estimated London population size.}
#' }
#'
#' @usage data(londonbills)
#' @examples
#' data(londonbills)
#' subset(londonbills, outbreak == 1665) # Great Plague of London
#'
#' @name londonbills
"londonbills"

#' Latent period distribution (pneumonic plague)
#'
#' @description
#' Distributions (observed and fitted) of the latent period
#' of pneumonic plague.
#'
#' @format
#' A data frame with 13 rows and 3 variables:
#'
#' \describe{
#'   \item{`days`}{\[integer\] Latent period in days, from 0 to 12.}
#'   \item{`observed`}{\[numeric\] Observed relative frequency
#'     out of 224 cases.
#'   }
#'   \item{`fitted`}{\[numeric\] Probability from a fitted lognormal
#'     distribution.
#'   }
#' }
#'
#' @source
#' Obtained from Figure 1A in \insertCite{GaniLeac04;textual}{epigrowthfit}.
#'
#' @references
#' \insertRef{GaniLeac04}{epigrowthfit}
#'
#' @usage data(plague_latent_period)
#' @name plague_latent_period
"plague_latent_period"

#' Infectious period distribution (pneumonic plague)
#'
#' @description
#' Distributions (observed and fitted) of the infectious period
#' of pneumonic plague.
#'
#' @format
#' A data frame with 13 rows and 3 variables:
#'
#' \describe{
#'   \item{`days`}{\[integer\] Infectious period in days, from 0 to 12.}
#'   \item{`observed`}{\[numeric\] Observed relative frequency
#'     out of 225 cases.
#'   }
#'   \item{`fitted`}{\[numeric\] Probability from a fitted lognormal
#'     distribution.
#'   }
#' }
#'
#' @source
#' Obtained from Figure 1B in \insertCite{GaniLeac04;textual}{epigrowthfit}.
#'
#' @references
#' \insertRef{GaniLeac04}{epigrowthfit}
#'
#' @usage data(plague_infectious_period)
#' @name plague_infectious_period
"plague_infectious_period"
