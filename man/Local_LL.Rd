\name{Local_LL}
\alias{Local_LL}
\title{Value of the Log-Likelihood Function L, where Input is in Phi-Parametrization}
\description{
  Gives the value of 
  
  \deqn{L(\varphi) = \sum_{i=1}^m w_i \varphi(x_i) - \int_{x_1}^{x_m} \exp(\varphi(t)) dt.}
}
\usage{Local_LL(x, w, phi)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}}}, where \eqn{w_1 > 0} and \eqn{w_m > 0}. These raw weights are 
            normalized in order to sum to one. Default: \eqn{w_i = 1/m}.}
  \item{phi}{Some vector \eqn{{\bold{\varphi}}} of the same length as \eqn{{\bold{x}}} and \eqn{{\bold{w}}}.}
}
\value{Value \eqn{L=L({\bold{\varphi}})} of the log-likelihood function is returned.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@freesurf.ch}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}

\keyword{htest}
\keyword{nonparametric}


