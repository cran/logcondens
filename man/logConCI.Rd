\name{logConCI}
\alias{logConCI}
\title{Compute pointwise confidence interval for a density assuming log-concavity}
\description{Compute approximate confidence interval for the true log-concave density, on a grid of points. 
The confidence interval at a fixed point is based on the pointwise asymptotic theory for the log-concave maximum likelihood estimator (MLE) developed in Balabdaoui, Rufibach, and Wellner (2009).}
\usage{logConCI(res, xx0, conf.level = c(0.8, 0.9, 0.95, 0.99)[3], 
type = c("DR", "ks", "nrd")[2], htype = c("hscv", "hlscv", "hpi", "hns")[4])}
\arguments{
    \item{res}{An object of class \code{dlc}, usually a result of a call to \code{logConDens}.}
    \item{xx0}{Vector of grid points at which to calculate the confidence interval.}
    \item{conf.level}{Confidence level for the confidence interval(s).}
    \item{type}{Vector of strings indicating estimating method for nuisance parameters.}	
    \item{htype}{Vector of strings indicating bandwidth selection method if \code{type = ks}.}
}
\value{The function returns a list containing the following elements:
    \item{fhat}{MLE evaluated at grid points.}
    \item{up_DR}{Upper confidence interval limit when \code{type = DR}.}
    \item{lo_DR}{Lower confidence interval limit when \code{type = DR}.}
    \item{up_ks_hscv}{Upper confidence interval limit when \code{type = ks} and \code{htype = hscv}.}
    \item{lo_ks_hscv}{Lower confidence interval limit when \code{type = ks} and \code{htype = hscv}.}
    \item{up_ks_hlscv}{Upper confidence interval limit when \code{type = ks} and \code{htype = hlscv}.}
    \item{lo_ks_hlscv}{Lower confidence interval limit when \code{type = ks} and \code{htype = hlscv}.}
    \item{up_ks_hpi}{Upper confidence interval limit when \code{type = ks} and \code{htype = hpi}.}
    \item{lo_ks_hpi}{Lower confidence interval limit when \code{type = ks} and \code{htype = hpi}.}
    \item{up_ks_hns}{Upper confidence interval limit when \code{type = ks} and \code{htype = hns}.}
    \item{lo_ks_hns}{Lower confidence interval limit when \code{type = ks} and \code{htype = hns}.}
    \item{up_nrd}{Upper confidence interval limit when \code{type = nrd}.}
    \item{lo_nrd}{Lower confidence interval limit when \code{type = nrd}.}
}
\details{In Balabdaoui, Rufibach, and Wellner (2009) it is shown that (if the true density is strictly log-concave) the limiting distribution of the MLE of a log-concave density \eqn{\widehat f_n} at a point \eqn{x} is

\deqn{n^{2/5}(\widehat f_n(x)-f(x)) \to c_2(x) \bar{C}(0).}

The nuisance parameter \eqn{c_2(x)} depends on the true density \eqn{f} and the second derivative of its logarithm. The limiting process \eqn{\bar{C}(0)} is found as the second derivative at zero of a particular operator (called the "envelope") of an integrated Brownian motion plus \eqn{t^4}.   

The confidence intervals are based on inverting the above limit using estimated quantiles of \eqn{\bar{C}(0)}, and estimating the nuisance parameter \eqn{c_2(x)}. The options for the function \code{logConCI} provide different ways to estimate this nuisance parameter. If \code{type="DR"}, \eqn{c_2(x)} is estimated using derivatives of the smoothed MLE as calculated by the function \code{logConDens} (this method does not perform well in simulations and is therefore not recommended). If \code{type="ks"}, \eqn{c_2(x)} is estimated using kernel density estimates of the true density and its first and second derivatives. This is done using the \code{R} package \pkg{ks}, and, with this option, a bandwidth selection method \code{htype} must also be chosen. The choices in \code{htype} correspond to the various options for bandwidth selection available in \pkg{ks}. If \code{type="nrd"}, the second derivative of the logarithm of the true density in \eqn{c_2(x)} is estimated assuming a normal reference distribution. The default option is \code{type="ks"} with \code{htype="hns"}.  

Currently available confidence levels are 80\%, 90\%, 95\% and 99\%, with a default of 95\%.
}
\references{
Azadbakhsh, M., Jankowski, H. and Gao, X. (2012)
Computing confidence intervals for log-concave densities. 
Preprint. 

Baladbaoui, F., Rufibach, K. and Wellner, J. (2009)
Limit distribution theory for maximum likelihood estimation of a log-concave density.
\emph{Ann. Statist.}, \bold{37(3)}, 1299--1331.

Tarn Duong (2012). ks: Kernel smoothing. 
R package version 1.8.10. \url{http://CRAN.R-project.org/package=ks}
}
\author{
Mahdis Azadbakhsh 

Hanna Jankowski, \email{hkj@mathstat.yorku.ca}, \cr \url{http://www.math.yorku.ca/~hkj/}
}
\examples{

## ===================================================
## Confidence intervals at a fixed point for the density
## ===================================================
data(reliability)
x.rel <- sort(reliability)

# calculate 95% confidence interval and plot the result
grid <- seq(min(x.rel), max(x.rel), length.out = 200)
res <- logConDens(x.rel)
ci  <- logConCI(res, grid, type = "nrd")	

par(las = 1, mar = c(2.5, 3.5, 0.5, 0.5))
hist(x.rel, n = 25, col = gray(0.9), main = "", freq = FALSE, 
    xlab = "", ylab = "", ylim = c(0, 0.0065), border = gray(0.5))
lines(grid, ci$fhat, col = "red", lwd = 2)
lines(grid, ci$lo_nrd, col = "red", lwd = 2, lty = 2)
lines(grid, ci$up_nrd, col = "red", lwd = 2, lty = 2)
}
\keyword{htest}
\keyword{nonparametric}
