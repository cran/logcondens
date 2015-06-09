\name{Local_LL}
\alias{Local_LL}
\title{Value of the Log-Likelihood Function L, where Input is in Phi-Parametrization}
\description{
  Gives the value of 
  
  \deqn{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - \int_{x_1}^{x_m} \exp(\phi(t)) dt.}{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - int_{x_1}^{x_m} \exp(\phi(t)) dt.}
}
\usage{Local_LL(x, w, phi)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}_m}}{x_m}.}
  \item{phi}{Some vector \eqn{{\bold{\phi}}}{\phi} of the same length as \eqn{{\bold{x}}}{x} and \eqn{{\bold{w}}}{w}.}
}
\value{Value \eqn{L=L({\bold{\phi}})}{L=L(\phi)} of the log-likelihood function is returned.}
\note{This function is not intended to be invoked by the end user.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/content/staff/personalhomepages/duembgen/index_eng.html}}

\keyword{htest}
\keyword{nonparametric}


