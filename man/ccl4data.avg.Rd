\name{ccl4data.avg}
\alias{ccl4data.avg}

\title{Closed chamber study of CCl4 metabolism by rats.}
\usage{data(ccl4data.avg)}
\description{
  The results of a closed chamber experiment to determine metabolic
  parameters for CCl4 (carbon tetrachloride) in rats.  This is the
  summary version of \code{\link{ccl4data}}.  Each record is the average
  for the time point of all animals exposed to the given initial chamber
  concentration.
  
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{time}{The time of the observation (hours after starting)}
    \item{initconc}{Initial chamber concentration (in ppm)}
    \item{ChamberConc}{Mean chamber concentration at the specified time (ppm)}
  }
}
\source{
 Evans, et al. 1994 Applications of sensitivity analysis to a
physiologically based pharmacokinetic model for carbon tetrachloride in
rats.  Toxicology and Applied Pharmacology 128: 36 -- 44.

}
\examples{
data(ccl4data.avg)
}
\keyword{datasets}
