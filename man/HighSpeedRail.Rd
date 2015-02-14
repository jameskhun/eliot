\name{HighSpeedRail}
\alias{HighSpeedRail}
\docType{data}
\title{
2010 and 2011 US Federal High Speed Rail Grant Data}
\description{
Data on high speed rail grants from the US federal government to the US states, in 2010 and 2011.}
\usage{data(HighSpeedRail)}
\format{
  A data frame with 50 observations on the following 5 variables.
  \describe{
    \item{\code{Grant2010}}{a numeric vector}
    \item{\code{Grant2011}}{a numeric vector}
    \item{\code{LogPopulation}}{a numeric vector}
    \item{\code{PresVote2008}}{a numeric vector}
    \item{\code{RepublicanGov2011}}{a logical vector}
  }
}
\details{
This is a small data set of high speed rail grants from the US federal government to the US states, in 2010 and 2011.  Several Republican state governors rejected federal high speed rail grants between 2010 and 2011.  Many states receive zero dollars of grants, making this data suitable for analysis with a compound poisson gamma regression model.  2010 rail grant data downloaded in June 2011 from \url{http://www.whitehouse.gov/files/documents/100128\_1400-HSRAwards-Summary\_FRA$\%$20Revisions.pdf}.  2011 rail grant data downloaded June 2011 from \url{http://www.fra.dot.gov/rpd/HSIPR/462.shtml}.}
\source{
Lauderdale, B.E. (2012) "Compound Poisson-Gamma Regression Models for Dollar Outcomes that are Sometimes Zero", Political Analysis (forthcoming).
}
\examples{

\dontrun{
data(HighSpeedRail)
Estimates2010 <- CPGregression(Grant2010~LogPopulation+PresVote2008+RepublicanGov2011,data=HighSpeedRail)
Estimates2011 <- CPGregression(Grant2011~LogPopulation+PresVote2008+RepublicanGov2011,data=HighSpeedRail)
print(Estimates2010)
print(Estimates2011)
}

}
\keyword{datasets}
