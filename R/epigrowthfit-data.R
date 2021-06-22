#' Data in the \pkg{epigrowthfit} package
#'
#' Epidemic time series to which growth rates can be fit,
#' and a few other data sets.
#'
#' @details
#' Below is a list of available data sets with links to their
#' documentation:
#' \describe{
#' \item{\code{\link{canadacovid}}}{
#'   Daily confirmations of COVID-19 in Canadian provinces
#'   and territories, from first confirmation to May 8, 2021.
#' }
#' \item{\code{\link{covid_generation_interval}}}{
#'   Gamma distribution of the COVID-19 generation interval
#'   fit to data from a cluster of 45 cases in Tianjin, China.
#' }
#' \item{\code{\link{husting}}}{
#'   Counts of wills probated in the Court of Husting
#'   during four plague epidemics in 14th century London.
#' }
#' \item{\code{\link{canterbury}}}{
#'   Counts of wills probated in the Prerogative Court of Canterbury
#'   during 24 plague epidemics in 16th and 17th century London.
#' }
#' \item{\code{\link{londonparishes}}}{
#'   Weekly counts of burials listed in extant parish registers
#'   during 24 plague epidemics in 16th and 17th century London.
#' }
#' \item{\code{\link{londonbills}}}{
#'   Weekly counts of plague deaths recorded in the London Bills of
#'   Mortality during 24 plague epidemics in 16th and 17th century
#'   London.
#' }
#' \item{\code{\link{plague_latent_period}}}{
#'   Empirical distribution of the latent period of pneumonic plague.
#' }
#' \item{\code{\link{plague_infectious_period}}}{
#'   Empirical distribution of the infectious period of pneumonic plague.
#' }
#' }
#'
#' @name epigrowthfit-data
#' @keywords internal
NULL

#' COVID-19 in Canadian provinces and territories
#'
#' Daily confirmations of COVID-19 in Canadian provinces
#' and territories, from first confirmation to May 8, 2021.
#'
#' @format
#' A \link[=data.frame]{data frame} with 5451 rows and 4 variables:
#' \describe{
#' \item{\code{province}}{
#'   A \link{factor}. Canadian province or territory (postal abbreviation).
#' }
#' \item{\code{Date}}{
#'   A \link{Date} vector. Daily within each level of \code{province},
#'   except prior to invasion, when less frequent reports were common.
#' }
#' \item{\code{cases_tot}}{
#'   An \link{integer} vector. \code{cases_tot[i]} is the number of cases
#'   confirmed up to the end of \code{Date[i]} in \code{province[i]}.
#' }
#' \item{\code{cases_new}}{
#'   An \link{integer} vector. \code{cases_new[i]} is the number of cases
#'   confirmed from the end of \code{Date[i-1]} to the end of \code{Date[i]}
#'   in \code{province[i]}.
#' }
#' }
#'
#' @source
#' Raw data were downloaded from Michael Li's public
#' \href{https://github.com/wzmli}{Github repository}.
#' Up-to-date data can be downloaded
#' \href{https://wzmli.github.io/COVID19-Canada/COVID19_Canada.csv}{here}.
#'
#' @usage data(canadacovid)
#' @examples
#' data(canadacovid)
#' subset(canadacovid, province == "ON", -province)
#'
#' @name canadacovid
"canadacovid"

#' Husting wills
#'
#' Counts of wills probated in the Court of Husting
#' during four plague epidemics in 14th century London.
#'
#' @format
#' A \link[=data.frame]{data frame} with 723 rows and 4 variables:
#' \describe{
#' \item{\code{Date}}{
#'   A \link{Date} vector. Spacing varies, as only nonzero counts of wills
#'   are included.
#' }
#' \item{\code{wills}}{
#'   An \link{integer} vector. \code{wills[i]} is the number of wills
#'   written on \code{Date[i]}.
#' }
#' \item{\code{outbreak}}{
#'   An ordered \link{factor}, \link{split}ting the time series by plague
#'   outbreak. \link[=levels]{Levels} indicate the years in which outbreaks
#'   took place: 1348, 1361, 1368, and 1375.
#' }
#' \item{\code{severity}}{
#'   An ordered \link{factor}, classifying outbreaks as \code{"minor"}
#'   or \code{"major"}. All 14th century outbreaks are classified as major,
#'   as the data are too sparse to distinguish between minor and major
#'   outbreaks, in contrast with wills probated in the
#'   \link[=canterbury]{Prerogative Court of Canterbury}.
#' }
#' }
#'
#' @source
#' These data were transcribed from Sharpe (1889).
#'
#' @references
#' Sharpe RR. Calendar of wills proved and enrolled in the Court of Husting,
#' London, A.D. 1258--A.D. 1688. London: JC Francis; 1889.
#'
#' @usage data(husting)
#' @examples
#' data(husting)
#' subset(husting, outbreak == 1375, -outbreak)
#' @name husting
"husting"

#' Canterbury wills
#'
#' Counts of wills probated in the Prerogative Court of Canterbury
#' during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A \link[=data.frame]{data frame} with 32681 rows and 4 variables:
#' \describe{
#' \item{\code{Date}}{
#'   A \link{Date} vector. Spacing varies, as only nonzero counts of wills
#'   are included.
#' }
#' \item{\code{wills}}{
#'   An \link{integer} vector. \code{wills[i]} is the number of wills written
#'   on \code{Date[i]}.
#' }
#' \item{\code{outbreak}}{
#'   An ordered \link{factor}, \code{split}ting the time series by plague
#'   outbreak. \link[=levels]{Levels} indicate the years in which outbreaks
#'   took place: 1563, 1578, ..., 1647, and 1665 (the Great Plague of London).
#' }
#' \item{\code{severity}}{
#'   An ordered \link{factor}, classifying outbreaks as \code{"minor"}
#'   or \code{"major"}. An outbreak is classified as major if and only if
#'   plague deaths per week per 1000 individuals exceeded 5 at least once.
#' }
#' }
#'
#' @source
#' These data were retrieved from the National Archives (UK) in 2018
#' using this
#' \href{https://www.nationalarchives.gov.uk/help-with-your-research/research-guides/wills-1384-1858/}{research guide}.
#'
#' @usage data(canterbury)
#' @examples
#' data(canterbury)
#' subset(canterbury, outbreak == 1665, -outbreak) # Great Plague of London
#' @name canterbury
"canterbury"

#' Burials in London parish registers
#'
#' Weekly counts of burials listed in extant parish registers
#' during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A \link[=data.frame]{data frame} with 14198 rows and 4 variables:
#' \describe{
#' \item{\code{Date}}{
#'   A \link{Date} vector. Weekly, except for 340 instances of 8-day spacing.
#' }
#' \item{\code{burials}}{
#'   An \link{integer} vector. \code{burials[i]} is the number of burials
#'   registered from the end of \code{Date[i-1]} to the end of \code{Date[i]}.
#' }
#' \item{\code{outbreak}}{
#'   An ordered \link{factor}, \link{split}ting the time series by plague
#'   outbreak. \link[=levels]{Levels} indicate the years in which outbreaks
#'   took place: 1563, 1578, ..., 1647, and 1665 (the Great Plague of London).
#' }
#' \item{\code{severity}}{
#'   An ordered \link{factor}, classifying outbreaks as \code{"minor"}
#'   or \code{"major"}. An outbreak is classified as major if and only if
#'   plague deaths per week per 1000 individuals exceeded 5 at least once.
#' }
#' }
#'
#' @source
#' These data were retrieved by Cummins et al. with permission from
#' \href{https://www.ancestry.com/}{Ancestry.com}.
#' They have been made publicly available
#' \href{http://www.neilcummins.com}{here}.
#'
#'
#' @references
#' Cummins N, Kelly M, \'{O} Gr\'{a}da C. Living standards and plague
#' in London, 1560--1665. Econ Hist Rev. 2016;69:3--34.
#'
#' @usage data(londonparishes)
#' @examples
#' data(londonparishes)
#' subset(londonparishes, outbreak == 1665, -outbreak) # Great Plague of London
#' @name londonparishes
"londonparishes"

#' Plague deaths in the London Bills of Mortality
#'
#' Weekly counts of plague deaths recorded in the London Bills of
#' Mortality during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A \link[=data.frame]{data frame} with 2041 rows and 6 variables:
#' \describe{
#' \item{\code{Date}}{
#'   A \link{Date} vector. Weekly within each level of `outbreak`,
#'   except for 78 instances of spacing not equal to 7 days.
#' }
#' \item{\code{deaths_all_causes}}{
#'   An \link{integer} vector. \code{deaths_all_causes[i]} is the number
#'   of deaths due to all causes from the end of \code{Date[i-1]}
#'   to the end of \code{Date[i]}.
#' }
#' \item{\code{deaths_plague}}{
#'   An integer vector. \code{deaths_plague[i]} is the number
#'   of deaths due to plague from the end of \code{Date[i-1]}
#'   to the end of \code{Date[i]}.
#' }
#' \item{\code{population}}{
#'   An \link{integer} vector. Estimated London population size.
#' }
#' \item{\code{outbreak}}{
#'   An ordered \link{factor}, \link{split}ting the time series by plague
#'   outbreak. \link[=levels]{Levels} indicate the years in which outbreaks
#'   took place: 1563, 1578, ..., 1647, and 1665 (the Great Plague of London).
#' }
#' \item{\code{severity}}{
#'   An ordered \link{factor}, classifying outbreaks as \code{"minor"}
#'   or \code{"major"}. An outbreak is classified as major if and only if
#'   plague deaths per week per 1000 individuals exceeded 5 at least once.
#' }
#' }
#'
#' @usage data(londonbills)
#' @examples
#' data(londonbills)
#' subset(londonbills, outbreak == 1665, -outbreak) # Great Plague of London
#'
#' @name londonbills
"londonbills"

#' Latent period distribution (pneumonic plague)
#'
#' Empirical distribution of the latent period of pneumonic plague.
#'
#' @format
#' A \link[=data.frame]{data frame} with 12 rows and 2 variables:
#' \describe{
#' \item{\code{days}}{Latent period in days, from 1 to 12.}
#' \item{\code{relfreq}}{Observed relative frequency out of 224 cases.}
#' }
#'
#' @source
#' Obtained from Figure 1A in Gani and Leach (2004).
#'
#' @references
#' Gani R, Leach S. Epidemiological determinants for modeling
#' pneumonic plague outbreaks. Emerg Infect Dis. 2004;10:608--14.
#'
#' @usage data(plague_latent_period)
#' @examples
#' data(plague_latent_period)
#' plot(relfreq ~ days, data = plague_latent_period)
#'
#' @name plague_latent_period
"plague_latent_period"

#' Infectious period distribution (pneumonic plague)
#'
#' Empirical distribution of the infectious period of pneumonic plague.
#'
#' @format
#' A \link[=data.frame]{data frame} with 12 rows and 2 variables:
#' \describe{
#' \item{\code{days}}{Infectious period in days, from 1 to 12.}
#' \item{\code{relfreq}}{Observed relative frequency out of 225 cases.}
#' }
#'
#' @source
#' Obtained from Figure 1B in Gani and Leach (2004).
#'
#' @references
#' Gani R, Leach S. Epidemiological determinants for modeling
#' pneumonic plague outbreaks. Emerg Infect Dis. 2004;10:608--14.
#'
#' @usage data(plague_infectious_period)
#' @examples
#' data(plague_infectious_period)
#' plot(relfreq ~ days, data = plague_infectious_period)
#'
#' @name plague_infectious_period
"plague_infectious_period"

#' Generation interval distribution (COVID-19)
#'
#' Gamma distribution of the COVID-19 generation interval fitted
#' to data from a cluster of 45 cases in Tianjin, China.
#'
#' @format
#' A \link{list} with 4 elements:
#' \describe{
#' \item{\code{shape}}{
#'   Estimated shape parameter. See \code{\link{dgamma}}.
#' }
#' \item{\code{scale}}{
#'   Estimated scale parameter. See \code{\link{dgamma}}.
#' }
#' \item{\code{breaks}}{
#'   An \link{integer} vector listing numbers of days. Equal to `0:20`.
#' }
#' \item{\code{probs}}{
#'   A \link[=double]{numeric} vector of length \code{\link{length}(breaks-1)}.
#'   \code{probs[i]} is the probability that the generation interval
#'   is between \code{breaks[i]} and \code{breaks[i+1]} days, conditional
#'   on \code{shape} and \code{scale}. Equal to
#'   \code{\link{diff}(\link{pgamma}(breaks, shape = shape, scale = scale))}.
#' }
#' }
#'
#' @source
#' \code{shape} and \code{scale} were computed from the mean and standard
#' deviation reported in Ganyani et al. (2020), Table 4, Scenario 2.
#'
#' @references
#' Ganyani T, Kremer C, Chen D, Torneri A, Faes C, Wallinga J, et al.
#' Estimating the generation interval for coronavirus disease (COVID-19)
#' based on symptom onset data, March 2020. Euro Surveill. 2020;25:2000257.
#'
#' @usage data(plague_infectious_period)
#' @examples
#' data(covid_generation_interval)
#' x <- 10^seq(-2, log10(20), length.out = 150L)
#' fx <- with(covid_generation_interval,
#'   dgamma(x, shape = shape, scale = scale)
#' )
#' plot(x, fx, type = "l", ylim = c(0, max(fx[-1L])))
#'
#' @name covid_generation_interval
"covid_generation_interval"
