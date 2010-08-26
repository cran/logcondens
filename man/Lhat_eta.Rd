\name{Lhat_eta}
\alias{Lhat_eta}
\title{Value of the Log-Likelihood Function L, where Input is in Eta-Parametrization}
\description{
  Gives the value of 
  
  \deqn{L(\varphi) = \sum_{i=1}^m w_i \varphi(x_i) - \int_{x_1}^{x_m} \exp(\varphi(t)) dt}
  
  where \eqn{\varphi} is parametrized via 
  
  \deqn{{\bold{\eta}}({\bold{\varphi}}) = \Bigl(\varphi_1, \Bigl(\eta_1 + \sum_{j=2}^i (x_i-x_{i-1})\eta_i\Bigr)_{i=2}^m\Bigr).} 
}
\usage{Lhat_eta(x, w, eta)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}_m}}.}
  \item{eta}{Some vector \eqn{{\bold{\eta}}} of the same length as \eqn{{\bold{x}}} and \eqn{{\bold{w}}}.}
}
\value{Value \eqn{L({\bold{\varphi}}) = L({\bold{\varphi}}({\bold{\eta}}))} of the log-likelihood function is returned.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}

\keyword{htest}
\keyword{nonparametric}

