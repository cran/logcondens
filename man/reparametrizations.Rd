\name{reparametrizations}
\alias{reparametrizations}
\alias{phieta}
\alias{etaphi}
\title{Changes Between Parametrizations}
\description{
  Given a vector \eqn{(\varphi_1, \ldots, \varphi_m)} representing the values of a piecewise linear concave function at
  \eqn{x_1, \ldots, x_m,} \code{\link{etaphi}} returns a column vector with the entries
   
  \deqn{{\bold{\eta}} = \Bigl(\varphi_1, \Bigl(\eta_1 + \sum_{j=2}^m (x_i-x_{i-1})\eta_i\Bigr)_{i=2}^m\Bigr). }
  
  The function \code{\link{phieta}} returns a vector with the entries 
  
  \deqn{{\bold{\varphi}} = \Bigl(\eta_1, \Bigl(\frac{\varphi_i-\varphi_{i-1}}{x_i-x_{i-1}}\Bigr)_{i=2}^m\Bigr).}
}
\usage{etaphi(x, eta)
phieta(x, phi)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{eta}{Column vector with entries \eqn{\eta_i = \eta(x_i).}}
  \item{phi}{Column vector with entries \eqn{\varphi_i = \varphi(x_i).}}
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@freesurf.ch}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}

\keyword{htest}
\keyword{nonparametric}
