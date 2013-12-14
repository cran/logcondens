\name{Jfunctions}
\alias{Jfunctions}
\alias{J00}
\alias{J10}
\alias{J11}
\alias{J20}
\title{Numerical Routine J and Some Derivatives}
\description{
  J00 represents the function \eqn{J(x, y, v),} where for real numbers \eqn{x, y} and \eqn{v \in [0, 1],}
  
    \deqn{J(x, y, v) = \int_0^v \exp((1-t)x + ty) d t = \frac{\exp(x + v(y - x)) - \exp(x)}{y - x}.}{J(x, y, v) = int_0^v exp((1 - t) x + t y) d t = (exp(x + v(y - x)) - exp(x))/(y - x).}
    
  The functions Jab give the respective derivatives \eqn{J_{ab}} for \eqn{v = 1}, i.e.
  
    \deqn{J_{ab}(x, y) = \frac{\partial^{a+b}}{\partial x^a \partial y^b} J(x, y).}{J_{ab}(x, y) = (partial ^ {a + b}) / (\partial x ^ a \partial y ^ b) J(x, y).}
    
  Specifically, 
  
    \deqn{J_{10}(x, y) = \frac{\exp(y) - \exp(x) - (y - x) \exp(x)}{(y - x)^2};}{J_{10}(x, y) = (exp(y) - exp(x) - (y - x) exp(x))/((y - x) ^ 2);}

    \deqn{J_{11}(x, y) = \frac{(y - x)(\exp(x) + \exp(y)) + 2 (\exp(y) - \exp(x))}{(y - x)^3};}{J_{11}(x, y) = ((y - x)(exp(x) + exp(y)) + 2 (exp(y) - exp(x)))/((y - x) ^ 3);}
  
    \deqn{J_{20}(x, y) = 2\frac{\exp(y) - \exp(x) - (y - x)\exp(x)-(y - x)^2 \exp(x)}{(y - x)^3}.}{J_{20}(x, y) = 2(exp(y) - exp(x) - (y - x) exp(x) - (y - x) ^ 2 exp(x)) / ((y - x) ^ 3).}
}
\usage{
J00(x, y, v)
J10(x, y)
J11(x, y)
J20(x, y)}
\arguments{
  \item{x}{Vector of length \eqn{d} with real entries.}
  \item{y}{Vector of length \eqn{d} with real entries.}
  \item{v}{Number in \eqn{[0, 1]^d}.}
}
\value{Value of the respective function.}
\note{Taylor approximations are used if \eqn{y-x} is small. We refer to Duembgen et al (2011, Section 6) for 
details.

These functions are not intended to be invoked by the end user.}
\references{
Duembgen, L, Huesler, A. and Rufibach, K. (2010)
Active set and EM algorithms for log-concave densities based on complete and censored data. 
Technical report 61, IMSV, Univ. of Bern, available at \url{http://arxiv.org/abs/0707.4643}.

Duembgen, L. and Rufibach, K. (2011)
logcondens: Computations Related to Univariate Log-Concave Density Estimation. 
\emph{Journal of Statistical Software}, \bold{39(6)}, 1--28. \url{http://www.jstatsoft.org/v39/i06}
}

\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}

\keyword{htest}
\keyword{nonparametric}
