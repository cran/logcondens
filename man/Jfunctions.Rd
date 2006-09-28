\name{Jfunctions}
\alias{Jfunctions}
\alias{J00}
\alias{J10}
\alias{J11}
\alias{J20}
\title{Numerical Routine J and Some Derivatives, to be used by activeSetLogCon}
\description{
  J00 represents the function \eqn{J(x, y),} where for real numbers \eqn{x, y,}
  
    \deqn{J(x, y) = \frac{\exp(y) - \exp(x)}{y - x}.}
    
  Jab gives the respective derivatives, i.e.
  
    \deqn{J10(x, y) = \frac{\exp(y) - \exp(x) - (y - x) \exp(x)}{(y - x)^2};}

    \deqn{J11(x, y) = \frac{(y - x)(\exp(x) + \exp(y)) + 2 (\exp(y) - \exp(x))}{(y - x)^3};}
  
    \deqn{J20(x, y) = 2\frac{\exp(y) - \exp(x) - (y - x)\exp(x)-(y - x)^2 \exp(x)}{(y - x)^3}.}
}
\usage{J00(x, y)
J10(x, y)
J11(x, y)
J20(x, y)}
\arguments{
  \item{x}{Arbitrary real number.}
  \item{y}{Arbitrary real number.}
}
\value{Value of the respective function.}
\note{Taylor approximations are used if \eqn{y-x} is small.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}

\keyword{htest}
\keyword{nonparametric}
