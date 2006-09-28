\name{logcon-package}
\alias{logcon-package}
\alias{logcon}
\docType{package}
\title{Estimate a Log-Concave Probability Density from i.i.d. Observations}
\description{Given a vector of observations \eqn{x = (x_1, \ldots, x_m)} with pairwise distinct entries and
a vector of weights \eqn{w=(w_1, \ldots, w_m)} s.t. \eqn{\sum_{i=1}^m w_i = 1}, this package allows to
 compute a concave, piecewise linear function \eqn{\widehat \varphi_m} on \eqn{[x_1, x_m]} with knots only in \eqn{\{x_1, \ldots, x_m\}}
 such that 
 
 \deqn{L(\varphi) = \sum_{i=1}^m w_i \varphi(x_i) - \int_{-\infty}^\infty \exp(\varphi(t)) dt}
 
 is maximal. According to the results in Duembgen and Rufibach (2006), this function \eqn{\widehat \varphi_m} maximizes the ordinary log-likelihood 
 
\deqn{\sum_{i=1}^m w_i \varphi(x_i)}
 
under the constraint that \eqn{\varphi} is concave. The corresponding function \eqn{\widehat f_m := \exp(\widehat\varphi_m)} is a log-concave probability density. Two algorithms are offered: An active set algorithm and one based on the pool-adjacent-violaters algorithm.}
\details{
\tabular{ll}{
Package: \tab logcondens\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-09-28\cr
License: \tab GPL version 2 or newer \cr
}
Use this package to estimate a log-concave density from i.i.d. observations. The package offers two algorithms which are implemented in the functions \code{\link{icmaLogCon}} and \code{\link{activeSetLogCon}}. Additionally, the value of the log-density, the density and the distribution function can be computed using \code{\link{evaluateLogConDens}}. To compute quantiles of the estimated distribution function use \code{\link{quantilesLogConDens}}.
}
\author{
Kaspar Rufibach (maintainer), \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}
}
\references{
Duembgen, L. and Rufibach, K. (2006) 
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
Preprint, IMSV, University of Bern.
              
Duembgen, L, Huesler, A. and Rufibach, K. (2006) 
Active set and EM algorithms for log-concave densities based on complete and censored data. 
Preprint, IMSV, University of Bern.

Rufibach, K. (2006a) 
Computing maximum likelihood estimators of a log-concave density function. 
To appear in \emph{Journal of Statistical Computation and Simulation.}
              
Rufibach K. (2006b) \emph{Log-concave Density Estimation and Bump Hunting for i.i.d. Observations.}
PhD Thesis, University of Bern, Switzerland and Georg-August University of Goettingen, Germany, 2006.
}

\examples{
## estimate gamma density
set.seed(1977)
x <- sort(rgamma(200, 2, 1))
res <- activeSetLogCon(x, w = NA, print = FALSE)

## plot resulting functions
par(mfrow = c(2, 2), mar = c(3, 2, 1, 2))
plot(x, exp(res$phi), type = 'l'); rug(x)
plot(x, res$phi, type = 'l'); rug(x)
plot(x, res$F, type = 'l'); rug(x)
plot(x, res$H, type = 'l'); rug(x)

## compute function values at an arbitrary point
x0 <- (x[100] + x[101]) / 2
evaluateLogConDens(x0, x, exp(res$phi), res$F, res$IsKnot)

## compute 0.5 quantile of estimated F
quantilesLogConDens(0.5, x, exp(res$phi), res$F, res$IsKnot)
}
\keyword{htest}
\keyword{nonparametric}
