\name{ccl4data}
\alias{ccl4data}
\non_function{}
\title{ Closed chamber study of CCl4 metabolism by rats. }
\usage{data(ccl4data)}
\description{
  The results of a closed chamber experiment to determine metabolic
  parameters for CCl4 (carbon tetrachloride) in rats.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{time}{The time (in hours after starting the experiment) }
    \item{initconc}{initial chamber concentration (ppm)}
    \item{animal}{This is a repeated measures design; this variable
      indicates which animal the observation pertains to}
    \item{ChamberConc}{chamber concentration at \code{time}, in ppm}
  }
}
\source{
Evans, et al. 1994 Applications of sensitivity analysis to a
physiologically based pharmacokinetic model for carbon tetrachloride in
rats.  Toxicology and Applied Pharmacology 128: 36 -- 44.
}
\examples{
data(ccl4data)
}
\keyword{datasets}
