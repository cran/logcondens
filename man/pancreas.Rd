\name{pancreas}
\alias{pancreas}
\docType{data}
\title{Data from pancreatic cancer serum biomarker study}
\description{First published by Wieand et al (1989), this dataset contains data on serum measurements for a cancer antigen (CA-125) and
a carbohydrate antigen (CA19.9) both of which are measured on a continuous non-negative scale. The measurements were taken within a case-control
study on 90 cases with pancreatic cancer and 51 controls who did not have cancer but pancreatitis. The primary question of the study was which one
of the two biomarkers best distinguishes cases from controls.}
\usage{data(pancreas)}
\format{
  A data frame with 141 observations on the following 3 variables.
  \describe{
    \item{\code{ca199}}{CA19.9 measurements.}
    \item{\code{ca125}}{CA125 measurements.}
    \item{\code{status}}{Patient status, 0 for controls and 1 for cases.}
  }
}
\source{The data was downloaded from http://labs.fhcrc.org/pepe/book/ on February 2, 2011.}
\references{
Wieand, S., Gail, M. H., James, B. R., and James, K.L. (1989).
A family of nonparametric statistics for comparing diagnostic markers with paired or unpaired data.
\emph{Biometrika}, \bold{76}, 585--592.

Pepe, M.S. (2003) \emph{The statistical evaluation of medical tests for classification and prediction}.
Oxford: Oxford University Press (Section 1.3.3).
}
\seealso{This data is analyzed in the help file for the function \code{\link{logConROC}}.}
\keyword{datasets}

