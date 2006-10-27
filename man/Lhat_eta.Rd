\name{Lhat_eta}
\alias{Lhat_eta}
\title{Value of the Log-Likelihood Function L, where Input is in Eta-Parametrization}
\description{
  Gives the value of 
  
  \deqn{L(\varphi) = \sum_{i=1}^m w_i \varphi(x_i) - \int_{x_1}^{x_m} \exp(\varphi(t)) dt}
  
  where \eqn{\varphi} is parametrized via 
  
  \deqn{\eta(\varphi) = \Bigl(\varphi_1, \Bigl(\eta_1 + \sum_{j=2}^i (x_i-x_{i-1})\eta_i\Bigr)_{i=2}^m\Bigr).} 
}
\usage{Lhat_eta(x, w, eta)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{x}, where \eqn{w_1 > 0} and \eqn{w_m > 0}. These raw weights are 
            normalized in order to sum to one. Default: \eqn{w_i = 1/m}.}
  \item{eta}{Some vector \eqn{\eta} of the same length as \eqn{x} and \eqn{w}.}
}
\value{Value \eqn{L(\varphi) = L(\varphi(\eta))} of the log-likelihood function is returned.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}

\keyword{htest}
\keyword{nonparametric}
