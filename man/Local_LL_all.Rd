\name{Local_LL_all}
\alias{Local_LL_all}
\title{Log-likelihood, New Candidate and Directional Derivative for L}
\description{
Computes the value of the log-likelihood function

  \deqn{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - \int_{x_1}^{x_m} \exp(\phi(t)) dt,}{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - int_{x_1}^{x_m} exp(\phi(t)) dt,}

a new candidate for \eqn{\phi} via the Newton method as well as the directional derivative of \eqn{{\bold{\phi}} \to L({\bold{\phi}})}{\phi \to L(\phi)} 
into that direction.
}
\usage{Local_LL_all(x, w, phi)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}_m}}{x_m}.}
  \item{phi}{Some vector \eqn{{\bold{\phi}}}{\phi} of the same length as \eqn{{\bold{x}}}{x} and \eqn{{\bold{w}}}{w}.}
}
\value{
  \item{ll}{Value \eqn{L(\phi)} of the log-likelihood function at \eqn{\phi.}}
  \item{phi_new}{New candidate for \eqn{\phi} via the Newton-method, using the complete Hessian matrix.}
  \item{dirderiv}{Directional derivative of \eqn{\phi \to L(\phi)} into the direction \eqn{\phi_{new}.}}
}
\note{This function is not intended to be invoked by the end user.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.imsv.unibe.ch/content/staff/personalhomepages/duembgen/index_eng.html}}

\keyword{htest}
\keyword{nonparametric}
