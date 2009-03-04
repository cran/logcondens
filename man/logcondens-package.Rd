\name{logcon-package}
\alias{logcon-package}
\alias{logcon}
\docType{package}
\title{Estimate a Log-Concave Probability Density from iid Observations}
\description{Given a vector of observations \eqn{{\bold{x}} = (x_1, \ldots, x_m)} with pairwise distinct entries and
a vector of weights \eqn{{\bold{w}}=(w_1, \ldots, w_m)} s.t. \eqn{\sum_{i=1}^m w_i = 1}, this package allows to compute a 
concave, piecewise linear function \eqn{\widehat \varphi_m} on \eqn{[x_1, x_m]} with knots only in 
\eqn{\{x_1, \ldots, x_m\}} such that 
 
 \deqn{L(\varphi) = \sum_{i=1}^m w_i \varphi(x_i) - \int_{-\infty}^\infty \exp(\varphi(t)) dt}
 
 is maximal. According to the results in Duembgen and Rufibach (2009), this function \eqn{\widehat \varphi_m} maximizes the ordinary log-likelihood 
 
\deqn{\sum_{i=1}^m w_i \varphi(x_i)}
 
under the constraint that \eqn{\varphi} is concave and that the corresponding density function \eqn{\varphi}
integrates to one. Two algorithms are offered: An active set and an iterative algorithm based on the pool-adjacent-violaters algorithm.}
\details{
\tabular{ll}{
Package: \tab logcondens\cr
Type: \tab Package\cr
Version: \tab 1.3.3 \cr
Date: \tab 2009-03-04\cr
License: \tab GPL version 2 or newer \cr
}
Use this package to estimate a log-concave density from i.i.d. observations. The package offers two algorithms which are 
implemented in the functions \code{\link{icmaLogCon}} and \code{\link{activeSetLogCon}}. Additionally, the value of the 
log-density, the density and the distribution function can be computed using \code{\link{evaluateLogConDens}}. 
To compute quantiles of the estimated distribution function use \code{\link{quantilesLogConDens}}.
}
\author{
Kaspar Rufibach (maintainer), \email{kaspar.rufibach@gmail.com} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}

Kaspar Rufibach gratefully acknowledges support by the Swiss National Science Foundation SNF, \url{http://www.snf.ch}.

MatLab code with an implementation of the active set algorithm is available on request from Lutz Duembgen.}
\references{
Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.  

Duembgen, L, Huesler, A. and Rufibach, K. (2007)
Active set and EM algorithms for log-concave densities based on complete and censored data.
Technical report 61, IMSV, Univ. of Bern, available at \url{http://arxiv.org/abs/0707.4643}.

Rufibach K. (2006) \emph{Log-concave Density Estimation and Bump Hunting for i.i.d. Observations.}
PhD Thesis, University of Bern, Switzerland and Georg-August University of Goettingen, Germany, 2006.
\cr Available at \url{http://www.stub.unibe.ch/download/eldiss/06rufibach_k.pdf}.

Rufibach, K. (2007)
Computing maximum likelihood estimators of a log-concave density function.
\emph{J. Stat. Comput. Simul.} \bold{77}, 561--574.
}

\examples{
## estimate gamma density
set.seed(1977)
x <- sort(rgamma(100, 2, 1))
res <- activeSetLogCon(x, w = NA, print = TRUE)

## compare performance to ICMA
res2 <- icmaLogCon(x, w = NA, T1 = 2000, robustif = TRUE, print = TRUE)

res$L
res2$L

## plot resulting functions
par(mfrow = c(2, 2), mar = c(3, 2, 1, 2))
plot(x, exp(res$phi), type = 'l'); rug(x)
plot(x, res$phi, type = 'l'); rug(x)
plot(x, res$Fhat, type = 'l'); rug(x)
plot(x, res$H, type = 'l'); rug(x)

## compute function values at an arbitrary point
x0 <- (x[50] + x[51]) / 2
evaluateLogConDens(x0, x, res$phi, res$Fhat, res$IsKnot)

## compute 0.5 quantile of Fhat
quantilesLogConDens(0.5, x, res$phi, res$Fhat)
}
\keyword{htest}
\keyword{nonparametric}
