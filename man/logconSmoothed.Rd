\name{logconSmoothed}
\Rdversion{1.1}
\alias{logconSmoothed}
\title{Compute convolution of log-concave density estimator with a normal kernel}
\description{Compute the convolution of the log-concave density estimator as output by \code{\link{activeSetLogCon}}
with a normal kernel. Since log-concavity is preserved under convolution and the normal density is log-concave, 
the resulting smoothed log-concave density estimator is still log-concave.}
\usage{logconSmoothed(x, w, phi, gam = NULL, xs = NULL, sig = NA)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}}}, where \eqn{w_1 > 0} and \eqn{w_m > 0}. These raw weights are 
            normalized in order to sum to one. Default: \eqn{w_i = 1 / m.}}
  \item{phi}{Column vector with entries \eqn{\widehat \varphi_m(x_i)} as output by \code{\link{activeSetLogCon}}
             or \code{\link{icmaLogCon}} (where the output of the latter is \eqn{f}, the density itself).}
  \item{gam}{The variance of the normal kernel. If equal to \code{NULL}, \code{gam} is chosen
             such that the variances of the sample \eqn{x_1, \ldots, x_m} and \eqn{\widehat f_n^*} coincide.}
  \item{xs}{Either provide a vector of support points where the smoothed estimator should be computed at, or
  leave as \code{NULL}. Then, a sufficiently width equidistant grid of points will be used.}
  \item{sig}{The variance of the initial, i.e. unweighted sample. This number can not be reconstructed
  from \code{x} and \code{w} if we have unequal weights and must therefore be given as input (but only if
  \code{gam} should be chosen as described above).}
}
\details{In Duembgen and Rufibach (2009, Section 3) the log-concave density estimator is illustrated using
real-life data. This function can be used to compute the smoothed log-concave density estimator briefly described
in that Section and added as a competitor to Figure 2.}
\value{
\item{f.smoothed}{Column vector with entries \eqn{\widehat f_m(x_i)}.}
\item{gam}{Variance of the normal kernel.}
\item{xs}{Support points at which the smoothed estimator was computed.}
}
\references{
Duembgen, L, Huesler, A. and Rufibach, K. (2007)
Active set and EM algorithms for log-concave densities based on complete and censored data. 
Technical report 61, IMSV, Univ. of Bern, available at \url{http://arxiv.org/abs/0707.4643}.

Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.  
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}
}
\examples{

## ===================================================
## Reproduce Fig. 2 in Duembgen & Rufibach (2009)
## ===================================================

## Set parameters
data(reliability)
x0 <- sort(reliability)
n0 <- length(x0)
x <- unique(x0)
w <- as.vector(table(x0))
w <- w / n0
n <- length(x)
res <- activeSetLogCon(x, w, print = TRUE)
phi <- res$phi
f <- exp(phi)

## compute normal density
mu <- mean(x0)
sig <- sd(x0)

## compute smoothed log-concave density estimator
f.smoothed <- logconSmoothed(x = x, w = w, phi, gam = NULL, sig = sig)
xs <- f.smoothed$xs

## compute kernel density
h <- sig / sqrt(n0)
f.kernel <- rep(NA, length(xs))
for (i in 1:length(xs)){
    xi <- xs[i]
    f.kernel[i] <- mean(dnorm(xi, mean = x0, sd = h))
}

## compute normal density
f.normal <- dnorm(xs, mean = mu, sd = sig)

## ===================================================
## Plot resulting densities, i.e. reproduce Fig. 2
## in Duembgen and Rufibach (2009)
## ===================================================
plot(0, 0, type = 'n', xlim = range(xs), ylim = c(0, 6.5 * 10^-3))
rug(x)
lines(x, f, col = 2)
lines(xs, f.normal, col = 3)
lines(xs, f.kernel, col = 4)
lines(xs, f.smoothed$f.smoothed, lwd = 3, col = 5)
legend("topleft", c("log-concave", "normal", "kernel", 
    "log-concave smoothed"), lty = 1, col = 2:5, bty = "n")


## ===================================================
## Plot log-densities
## ===================================================
plot(0, 0, type = 'n', xlim = range(xs), ylim = c(-20, -5))
legend("bottomright", c("log-concave", "normal", "kernel", 
    "log-concave smoothed"), lty = 1, col = 2:5, bty = "n")
rug(x)
lines(x, phi, col = 2)
lines(xs, log(f.normal), col = 3)
lines(xs, log(f.kernel), col = 4)
lines(xs, log(f.smoothed$f.smoothed), lwd = 3, col = 5)
}
\keyword{htest}
\keyword{nonparametric}
