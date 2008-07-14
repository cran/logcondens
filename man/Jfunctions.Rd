\name{Jfunctions}
\alias{Jfunctions}
\alias{J00}
\alias{J10}
\alias{J11}
\alias{J20}
\title{Numerical Routine J and Some Derivatives}
\description{
  J00 represents the function \eqn{J(x, y, v),} where for real numbers \eqn{x, y} and \eqn{v \in [0, 1],}
  
    \deqn{J(x, y, v) = \int_0^v \exp((1-t)x + ty) d t = \frac{\exp(x + v(y - x)) - \exp(x)}{y - x}.}
    
  The functions Jab give the respective derivatives \eqn{J_{ab}} for \eqn{v = 1}, i.e.
  
    \deqn{J_{ab}(x, y) = \frac{\partial^{a+b}}{\partial x^a \partial y^b} J(x, y).}
    
  Specifically, 
  
    \deqn{J_{10}(x, y) = \frac{\exp(y) - \exp(x) - (y - x) \exp(x)}{(y - x)^2};}

    \deqn{J_{11}(x, y) = \frac{(y - x)(\exp(x) + \exp(y)) + 2 (\exp(y) - \exp(x))}{(y - x)^3};}
  
    \deqn{J_{20}(x, y) = 2\frac{\exp(y) - \exp(x) - (y - x)\exp(x)-(y - x)^2 \exp(x)}{(y - x)^3}.}
}
\usage{J00(x, y, v)
J10(x, y)
J11(x, y)
J20(x, y)}
\arguments{
  \item{x}{Vector of length \eqn{d} with real entries.}
  \item{y}{Vector of length \eqn{d} with real entries.}
  \item{v}{Number in \eqn{[0, 1]^d}.}
}
\value{Value of the respective function.}
\note{Taylor approximations are used if \eqn{y-x} is small.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}

\keyword{htest}
\keyword{nonparametric}
