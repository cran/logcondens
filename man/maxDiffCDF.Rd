\name{maxDiffCDF}
\alias{maxDiffCDF}
\title{Compute maximal difference between CDFs corresponding to log-concave estimates}
\description{Compute the maximal difference between two estimated log-concave distribution functions, either
the MLEs or the smoothed versions. This function is used to set up a two-sample permutation test
that assesses the null hypothesis of equality of distribution functions.}
\usage{maxDiffCDF(res1, res2, which = c("MLE", "smooth"), n.grid = 500)}
\arguments{
  \item{res1}{An object of class \code{"dlc"}, usually a result of a call to \code{logConDens} for the first sample.}
  \item{res2}{An object of class \code{"dlc"}, usually a result of a call to \code{logConDens} for the second sample.}
  \item{which}{Indicate for which type of estimate the maximal difference should be computed.}
  \item{n.grid}{Number of grid points used to find zeros of \eqn{\hat f_{n_1}^* - \hat f_{n_2}^*} for the smooth estimate.}
}
\details{Given two i.i.d. samples \eqn{x_1, \ldots, x_{n_1}} and \eqn{y_1, \ldots, y_{n_2}} this function computes the
maxima of 

\deqn{D_1(t) = \hat F_{n_1}(t) - \hat F_{n_2}(t)}

and

\deqn{D_2(t) = \hat F^*_{n_1}(t) - \hat F^*_{n_2}(t).}
}
\value{
\item{test.stat}{A two-dimensional vector containing the above maxima.}
\item{location}{A two-dimensional vector where the maxima occur.}
}
\section{Warning}{Note that the algorithm that finds the maximal difference for the smoothed estimate is of approximative
nature only. It may fail for very large sample sizes.}
\references{
Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.  

Duembgen, L. and Rufibach, K. (2011)
logcondens: Computations Related to Univariate Log-Concave Density Estimation. 
\emph{Journal of Statistical Software}, \bold{39(6)}, 1--28. \url{http://www.jstatsoft.org/v39/i06}
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}
}
\examples{
n1 <- 100
n2 <- 120
x <- sort(rgamma(n1, 2, 1))
y <- sort(rgamma(n2, 2, 1))
res1 <- logConDens(x, smoothed = TRUE)
res2 <- logConDens(y, smoothed = TRUE)
d <- maxDiffCDF(res1, res2, n.grid = 200)

## log-concave estimate
xs <- seq(min(res1$xs, res2$xs), max(res1$xs, res2$xs), length = 200)
F1 <- matrix(NA, nrow = length(xs), ncol = 1); F2 <- F1
for (i in 1:length(xs)){
    F1[i] <- evaluateLogConDens(xs[i], res1, which = 3)[, "CDF"]
    F2[i] <- evaluateLogConDens(xs[i], res2, which = 3)[, "CDF"]
    }
par(mfrow = c(1, 2))
plot(xs, abs(F1 - F2), type = "l")
abline(v = d$location[1], lty = 2, col = 3)
abline(h = d$test.stat[1], lty = 2, col = 3)

## smooth estimate
xs <- seq(min(res1$xs, res2$xs), max(res1$xs, res2$xs), length = 200)
F1smooth <- matrix(NA, nrow = length(xs), ncol = 2); F2smooth <- F1smooth
for (i in 1:length(xs)){
    F1smooth[i, ] <- evaluateLogConDens(xs[i], res1, which = 4:5)[, 
        c("smooth.density", "smooth.CDF")]
    F2smooth[i, ] <- evaluateLogConDens(xs[i], res2, which = 4:5)[, 
        c("smooth.density", "smooth.CDF")]
    }
plot(xs, abs(F1smooth[, 2] - F2smooth[, 2]), type = "l")
abline(h = 0)
abline(v = d$location[2], lty = 2, col = c(3, 4))
abline(h = d$test.stat[2], lty = 2, col = c(3, 4))
}
\keyword{htest}
\keyword{nonparametric}
