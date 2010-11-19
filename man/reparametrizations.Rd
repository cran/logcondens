\name{reparametrizations}
\alias{reparametrizations}
\alias{phieta}
\alias{etaphi}
\title{Changes Between Parametrizations}
\description{
  Given a vector \eqn{(\phi_1, \ldots, \phi_m)} representing the values of a piecewise linear concave function at
  \eqn{x_1, \ldots, x_m,} \code{\link{etaphi}} returns a column vector with the entries
   
  \deqn{{\bold{\eta}} = \Bigl(\phi_1, \Bigl(\eta_1 + \sum_{j=2}^m (x_i-x_{i-1})\eta_i\Bigr)_{i=2}^m\Bigr).}{\eta = (\phi_1, (\eta_1 + \sum_{j=2}^m (x_i-x_{i-1})\eta_i)_{i=2}^m).}
  
  The function \code{\link{phieta}} returns a vector with the entries 
  
  \deqn{{\bold{\phi}} = \Bigl(\eta_1, \Bigl(\frac{\phi_i-\phi_{i-1}}{x_i-x_{i-1}}\Bigr)_{i=2}^m\Bigr).}{\phi = (\eta_1, ((\phi_i-\phi_{i-1})/(x_i-x_{i-1}))_{i=2}^m).}
}
\usage{etaphi(x, eta)
phieta(x, phi)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{eta}{Vector with entries \eqn{\eta_i = \eta(x_i).}}
  \item{phi}{Vector with entries \eqn{\phi_i = \phi(x_i).}}
}
\note{These functions are not intended to be invoked by the end user.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}

\keyword{htest}
\keyword{nonparametric}
