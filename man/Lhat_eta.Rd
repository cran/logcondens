\name{Lhat_eta}
\alias{Lhat_eta}
\title{Value of the Log-Likelihood Function L, where Input is in Eta-Parametrization}
\description{
  Gives the value of 
  
  \deqn{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - \int_{x_1}^{x_m} \exp(\phi(t)) dt}{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - int_{x_1}^{x_m} exp(\phi(t)) dt}
  
  where \eqn{\phi} is parametrized via 
  
  \deqn{{\bold{\eta}}({\bold{\phi}}) = \Bigl(\phi_1, \Bigl(\eta_1 + \sum_{j=2}^i (x_i-x_{i-1})\eta_i\Bigr)_{i=2}^m\Bigr).}{\eta(\phi) = (\phi_1, (\eta_1 + \sum_{j=2}^i (x_i-x_{i-1})\eta_i)_{i=2}^m).} 
}
\usage{Lhat_eta(x, w, eta)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}_m}}{x_m}.}
  \item{eta}{Some vector \eqn{{\bold{\eta}}}{\eta} of the same length as \eqn{{\bold{x}}}{x} and \eqn{{\bold{w}}}{w}.}
}
\value{Value \eqn{L({\bold{\phi}}) = L({\bold{\phi}}({\bold{\eta}}))}{L(\phi) = L(\phi(\eta))} of the log-likelihood function is returned.}
\note{This function is not intended to be invoked by the end user.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.imsv.unibe.ch/about_us/staff/prof_dr_duembgen_lutz/index_eng.html}}

\keyword{htest}
\keyword{nonparametric}

