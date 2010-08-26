\name{Q00}
\alias{Q00}
\title{Numerical Routine Q}
\description{This function is used in the computation of \eqn{\widehat f^*} and \eqn{\widehat F^*}.}
\usage{
Q00(x, a, u, v, gamma, QFhat = FALSE)}
\arguments{
  \item{x}{Number at which to compute \eqn{q} and/or \eqn{Q}.}
  \item{a}{Vector of length \eqn{m} with real entries.}
  \item{u}{Vector of length \eqn{m} with real entries.}
  \item{v}{Vector of length \eqn{m} with real entries.}
  \item{gamma}{Real number. Standard deviation to be used.}
  \item{QFhat}{Logical. Should \eqn{Q} be computed?}
}
\value{The vector(s) \eqn{q} and/or \eqn{Q}.}
\note{Taylor approximation is used if \eqn{a} is small.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}

\keyword{htest}
\keyword{nonparametric}
