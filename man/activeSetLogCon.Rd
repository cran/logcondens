\name{activeSetLogCon}
\alias{activeSetLogCon}
\alias{activeSet}
\title{Computes a Log-Concave Probability Density Estimate via an Active Set Algorithm}
\description{Given a vector of observations \eqn{{\bold{x}_n} = (x_1, \ldots, x_n)}{x_n = (x_1, \ldots, x_n)} 
with not necessarily equal entries, 
\code{\link{activeSetLogCon}} first computes vectors \eqn{{\bold{x}_m} = (x_1, \ldots, x_m)}{x_m = (x_1, \ldots, x_m)} 
and \eqn{{\bold{w}} = (w_1, \ldots, w_m)}{w = (w_1, \ldots, w_m)} where \eqn{w_i} is the weight of each \eqn{x_i} s.t. 
\eqn{\sum_{i=1}^m w_i = 1}.
Then, \code{\link{activeSetLogCon}} computes a concave, piecewise 
linear function \eqn{\widehat \phi_m}{\widehat \phi_m} on \eqn{[x_1, x_m]} with knots only in \eqn{\{x_1, \ldots, x_m\}}{{x_1, \ldots, x_m}} such that 
  
  \deqn{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - \int_{-\infty}^\infty \exp(\phi(t)) dt}{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - int_{-\infty}^\infty exp(\phi(t)) dt}
  
  is maximal. To accomplish this, an active set algorithm is used.} 
\usage{activeSetLogCon(x, xgrid = NULL, print = FALSE, w = NA)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, not necessarily unique.}
  \item{xgrid}{Governs the generation of weights for observations. See \code{\link{preProcess}} for details.}
  \item{print}{\code{print = TRUE} outputs the log-likelihood in every loop, \code{print = FALSE} does not. Make sure to tell \code{R} to output (press CTRL+W).}
  \item{w}{Optional vector of weights. If weights are provided, i.e. if \code{w != NA}, then \code{xgrid} is ignored.}
  }
\value{
  \item{xn}{Vector with initial observations \eqn{x_1, \ldots, x_n}.}  
  \item{x}{Vector of observations \eqn{x_1, \ldots, x_m} that was used to estimate the density. }
  \item{w}{The vector of weights that had been used. Depends on the chosen setting for \code{xgrid}.}
  \item{phi}{Vector with entries \eqn{\widehat \phi_m(x_i)}.}
  \item{IsKnot}{Vector with entries IsKnot\eqn{_i = 1\{\widehat \phi_m} has a kink at \eqn{x_i\}}.}
  \item{L}{The value \eqn{L(\widehat {\bold{\phi}}_m)}{L(\phi_m)} of the log-likelihood-function \eqn{L} at the 
    maximum \eqn{\widehat {\bold{\phi}}_m}{\widehat \phi_m}.}
  \item{Fhat}{A vector \eqn{(\widehat F_{m,i})_{i=1}^m} of the same size as \eqn{{\bold{x}}}{x} with entries 

  \deqn{\widehat F_{m,i} = \int_{x_1}^{x_i} \exp(\widehat \phi_m(t)) dt.}}
  \item{H}{Vector \eqn{(H_1, \ldots, H_m)'} where \eqn{H_i} is the derivative of 
  
  \deqn{t \to L(\phi + t\Delta_i)}
  
  at zero and \eqn{\Delta_i(x) = \min(x - x_i, 0).}}
  \item{n}{Number of initial observations.}  
  \item{m}{Number of unique observations.}  
  \item{knots}{Observations that correspond to the knots.}
  \item{mode}{Mode of the estimated density \eqn{\hat f_m}.}
  \item{sig}{The standard deviation of the initial sample \eqn{x_1, \ldots, x_n}.}
}
\references{
Duembgen, L, Huesler, A. and Rufibach, K. (2010)
Active set and EM algorithms for log-concave densities based on complete and censored data. 
Technical report 61, IMSV, Univ. of Bern, available at \url{http://arxiv.org/abs/0707.4643}.

Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68. 

Duembgen, L. and Rufibach, K. (2011)
logcondens: Computations Related to Univariate Log-Concave Density Estimation. 
\emph{Journal of Statistical Software}, \bold{39(6)}, 1--28. \url{http://www.jstatsoft.org/v39/i06}
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/content/staff/personalhomepages/duembgen/index_eng.html}}
\seealso{\code{\link{activeSetLogCon}} can be used to estimate a log-concave density. However, to generate an object of
class \code{dlc} that allows application of \code{\link{summary}} and \code{\link{plot}} we recommend to use \code{\link{logConDens}}.

The following functions are used by \code{\link{activeSetLogCon}}: 

\code{\link{J00}}, \code{\link{J10}}, \code{\link{J11}}, \code{\link{J20}},
\code{\link{Local_LL}}, \code{\link{Local_LL_all}}, \code{\link{LocalCoarsen}}, 
\code{\link{LocalConvexity}}, \code{\link{LocalExtend}}, \code{\link{LocalF}}, \code{\link{LocalMLE}}, 
\code{\link{LocalNormalize}}, \code{\link{MLE}}

Log concave density estimation via an iterative convex minorant algorithm can be performed using
\code{\link{icmaLogCon}}.
}
\examples{
## estimate gamma density
set.seed(1977)
n <- 200
x <- rgamma(n, 2, 1)
res <- activeSetLogCon(x, w = rep(1 / n, n), print = FALSE)

## plot resulting functions
par(mfrow = c(2, 2), mar = c(3, 2, 1, 2))
plot(res$x, exp(res$phi), type = 'l'); rug(x)
plot(res$x, res$phi, type = 'l'); rug(x)
plot(res$x, res$Fhat, type = 'l'); rug(x)
plot(res$x, res$H, type = 'l'); rug(x)

## compute and plot function values at an arbitrary point
x0 <- (res$x[100] + res$x[101]) / 2
Fx0 <- evaluateLogConDens(x0, res, which = 3)[, "CDF"]
plot(res$x, res$Fhat, type = 'l'); rug(res$x)
abline(v = x0, lty = 3); abline(h = Fx0, lty = 3)

## compute and plot 0.9-quantile of Fhat
q <- quantilesLogConDens(0.9, res)[2]
plot(res$x, res$Fhat, type = 'l'); rug(res$x)
abline(h = 0.9, lty = 3); abline(v = q, lty = 3)
}
\keyword{htest}
\keyword{nonparametric}
