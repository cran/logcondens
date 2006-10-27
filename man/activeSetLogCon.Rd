\name{activeSetLogCon}
\alias{activeSetLogCon}
\alias{activeSet}
\alias{log-concave}
\alias{logcondens}
\title{Computes a Log-Concave Probability Density Estimate via an Active Set Algorithm}
\description{Given a vector of observations \eqn{x = (x_1, \ldots, x_m)} with pairwise distinct entries and
a vector of weights \eqn{w=(w_1, \ldots, w_m)} s.t. \eqn{\sum_{i=1}^m w_i = 1}, \code{\link{activeSetLogCon}} computes a concave, piecewise 
linear function \eqn{\widehat \varphi_m} on \eqn{[x_1, x_m]} with knots only in \eqn{\{x_1, \ldots, x_m\}} such that 
  
  \deqn{L(\varphi) = \sum_{i=1}^m w_i \varphi(x_i) - \int_{-\infty}^\infty \exp(\varphi(t)) dt}
  
  is maximal. To accomplish this, an active set algorithm is used.} 
\usage{activeSetLogCon(x, w = NA, print = FALSE)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{x}, where \eqn{w_1 > 0} and \eqn{w_m > 0}. These raw weights are 
            normalized in order to sum to one. Default: \eqn{w_i = 1 / m.}}
  \item{print}{print = TRUE outputs log-likelihood in every loop, print = FALSE does not. Make sure to tell R to output (press CTRL+W).}
}
\value{
  \item{x}{Column vector with observations \eqn{x_1, \ldots, x_m}.}
  \item{phi}{Column vector with entries \eqn{\widehat \varphi_m(x_i)}.}
  \item{IsKnot}{Column vector with entries \cr
  IsKnot\eqn{_i = 1\{\widehat \varphi_m} has a kink at \eqn{x_i\}}.}
  \item{L}{The value \eqn{L(\widehat \varphi_m)} of the log-likelihood-function \eqn{L} at the maximum \eqn{\widehat \varphi_m}.}
  \item{F}{A vector \eqn{(\widehat F_{m,i})_{i=1}^m} of the same size as \eqn{x} with entries 
  
  \deqn{\widehat F_{m,i} = \int_{x_1}^{x_i} \exp(\widehat \varphi_m(t)) dt.}}
  
  \item{H}{Column vector \eqn{(H_1, \ldots, H_m)'} where \eqn{H_i} is the derivative of 
  
  \deqn{t \to L(\varphi + t\Delta)}
  
  at zero and \eqn{\Delta(x) = \min(x - x_k, 0).}}
}
\references{
Duembgen, L. and Rufibach, K. (2006)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
Preprint, IMSV, University of Bern.
              
Duembgen, L, Huesler, A. and Rufibach, K. (2006)
Active set and EM algorithms for log-concave densities based on complete and censored data. 
Preprint, IMSV, University of Bern.
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}
\seealso{The following functions are used by \code{\link{activeSetLogCon}}: 

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
x <- sort(rgamma(200, 2, 1))
res <- activeSetLogCon(x, w = NA, print = FALSE)

## plot resulting functions
par(mfrow = c(2, 2), mar = c(3, 2, 1, 2))
plot(x, exp(res$phi), type = 'l'); rug(x)
plot(x, res$phi, type = 'l'); rug(x)
plot(x, res$Fhat, type = 'l'); rug(x)
plot(x, res$H, type = 'l'); rug(x)

## compute and plot function values at an arbitrary point
x0 <- (x[100] + x[101]) / 2
Fx0 <- evaluateLogConDens(x0, x, res$phi, res$Fhat, res$IsKnot)[3]
plot(x, res$Fhat, type = 'l'); rug(x)
abline(v = x0, lty = 3); abline(h = Fx0, lty = 3)

## compute and plot 0.9-quantile of Fhat\coprod
q <- quantilesLogConDens(0.9, x, res$phi, res$Fhat)
plot(x, res$Fhat, type = 'l'); rug(x)
abline(h = 0.9, lty = 3); abline(v = q, lty = 3)
}
\keyword{htest}
\keyword{nonparametric}
