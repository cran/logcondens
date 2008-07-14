\name{Local_LL_all}
\alias{Local_LL_all}
\title{Log-likelihood, New Candidate and Directional Derivative for L}
\description{
Computes the value of the log-likelihood function

  \deqn{L(\varphi) = \sum_{i=1}^m w_i \varphi(x_i) - \int_{x_1}^{x_m} \exp(\varphi(t)) dt,}

a new candidate for \eqn{\varphi} via the Newton method as well as the directional derivative of \eqn{{\bold{\varphi}} \to L({\bold{\varphi}})} into that direction.
}
\usage{Local_LL_all(x, w, phi)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}}}, where \eqn{w_1 > 0} and \eqn{w_m > 0}. These raw weights are 
            normalized in order to sum to one. Default: \eqn{w_i = 1/m}.}
  \item{phi}{Some vector \eqn{{\bold{\varphi}}} of the same length as \eqn{{\bold{x}}} and \eqn{{\bold{w}}}.}
}
\value{
  \item{ll}{Value \eqn{L(\varphi)} of the log-likelihood function at \eqn{\varphi.}}
  \item{phi_new}{New candidate for \eqn{\varphi} via the Newton-method, using the complete Hessian matrix.}
  \item{dirderiv}{Directional derivative of \eqn{\varphi \to L(\varphi)} into the direction \eqn{\varphi_{new}.}}
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}

\keyword{htest}
\keyword{nonparametric}
