\name{logConCI}
\alias{logConCI}
\title{Compute pointwise confidence interval for a density assuming log-concavity}
\description{Compute approximate confidence interval for the true log-concave density, on a grid of points.  Two main approaches are implemented:  
In the first, the confidence interval at a fixed point is based on the pointwise asymptotic theory for the log-concave maximum likelihood estimator (MLE) developed in Balabdaoui, Rufibach, and Wellner (2009). In the second, the confidence interval is estimated via the boostrap.}
\usage{logConCI(res, xx0, conf.level = c(0.8, 0.9, 0.95, 0.99)[3], 
    type = c("DR", "ks", "nrd", "ECDFboot", "NPMLboot")[2], 
    htype = c("hscv", "hlscv", "hpi", "hns")[4], BB = 500)}
\arguments{
    \item{res}{An object of class \code{dlc}, usually a result of a call to \code{logConDens}.}
    \item{xx0}{Vector of grid points at which to calculate the confidence interval.}
    \item{conf.level}{Confidence level for the confidence interval(s). The default is 95\%.}
    \item{type}{Vector of strings indicating type of confidence interval to compute.  When \code{type = ks} is chosen, then \code{htype} should also be specified. The default is \code{type = ks}.}	
    \item{htype}{Vector of strings indicating bandwidth selection method if \code{type = ks}. The default is \code{htype = hns}.}
    \item{BB}{number of iterations in the bootstrap if \code{type = NPMLboot} or \code{type = ECDFboot}. The default is \code{BB = 500}.}
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
    \item{up_npml}{Upper confidence interval limit when \code{type = NPMLboot}.}
    \item{lo_npml}{Lower confidence interval limit when \code{boot = NPMLboot}.}
    \item{up_ecdf}{Upper confidence interval limit when \code{boot = ECDFboot}.}
    \item{lo_ecdf}{Lower confidence interval limit when \code{boot = ECDFboot}.}
}
\details{In Balabdaoui et al. (2009) it is shown that (if the true density is strictly log-concave) the limiting distribution of the MLE of a log-concave 
density \eqn{\widehat f_n} at a point \eqn{x} is

\deqn{n^{2/5}(\widehat f_n(x)-f(x)) \to c_2(x) \bar{C}(0).}

The nuisance parameter \eqn{c_2(x)} depends on the true density \eqn{f} and the second derivative of its logarithm. The limiting process \eqn{\bar{C}(0)} 
is found as the second derivative at zero of a particular operator (called the "envelope") of an integrated Brownian motion plus \eqn{t^4}.   

Three of the confidence intervals are based on inverting the above limit using estimated quantiles of \eqn{\bar{C}(0)}, and estimating the nuisance 
parameter \eqn{c_2(x)}. The options for the function \code{logConCI} provide different ways to estimate this nuisance parameter. If \code{type = "DR"}, 
\eqn{c_2(x)} is estimated using derivatives of the smoothed MLE as calculated by the function \code{logConDens} (this method does not perform well in 
simulations and is therefore not recommended). If \code{type="ks"}, \eqn{c_2(x)} is estimated using kernel density estimates of the true density and its 
first and second derivatives. This is done using the \code{R} package \pkg{ks}, and, with this option, a bandwidth selection method \code{htype} must also 
be chosen. The choices in \code{htype} correspond to the various options for bandwidth selection available in \pkg{ks}. If \code{type = "nrd"}, the second 
derivative of the logarithm of the true density in \eqn{c_2(x)} is estimated assuming a normal reference distribution. 

Two of the confidence intervals are based on the bootstrap. For \code{type = "ECDFboot"} confidence intervals based on re-sampling from the empirical 
cumulative distribution function are computed. For \code{type = "NPMLboot"} confidence intervals based on re-sampling from the nonparametric maximum 
likelihood estimate of log-concave density are computed. Bootstrap confidence intervals take a few minutes to compute!  

The default option is \code{type = "ks"} with \code{htype = "hns"}.  Currently available confidence levels are 80\%, 90\%, 95\% and 99\%, with a default 
of 95\%.   

Azadbakhsh et al. (2014) provides an empirical study of the relative performance of the various approaches available in this function.    
}
\references{
Azadbakhsh, M., Jankowski, H. and Gao, X. (2014).
Computing confidence intervals for log-concave densities. 
\emph{Comput. Statist. Data Anal.}, \bold{75}, 248--264.

Baladbaoui, F., Rufibach, K. and Wellner, J. (2009)
Limit distribution theory for maximum likelihood estimation of a log-concave density.
\emph{Ann. Statist.}, \bold{37(3)}, 1299--1331.

Tarn Duong (2012). ks: Kernel smoothing. 
R package version 1.8.10. \url{https://CRAN.R-project.org/package=ks}
}
\author{
Mahdis Azadbakhsh 

Hanna Jankowski, \email{hkj@mathstat.yorku.ca}, \cr \url{http://www.math.yorku.ca/~hkj/}
}
\examples{

\dontrun{
## ===================================================
## Confidence intervals at a fixed point for the density
## ===================================================
data(reliability)
x.rel <- sort(reliability)

# calculate 95% confidence interval(s) and plot the result
grid <- seq(min(x.rel), max(x.rel), length.out = 200)
res <- logConDens(x.rel)
ci  <- logConCI(res, grid, type = c("nrd", "ECDFboot"))	

par(las = 1, mar = c(2.5, 3.5, 0.5, 0.5))
hist(x.rel, n = 25, col = gray(0.9), main = "", freq = FALSE, 
    xlab = "", ylab = "", ylim = c(0, 0.0065), border = gray(0.5))
lines(grid, ci$fhat, col = "black", lwd = 2)
lines(grid, ci$lo_nrd, col = "red", lwd = 2, lty = 2)
lines(grid, ci$up_nrd, col = "red", lwd = 2, lty = 2)
lines(grid, ci$lo_ecdf, col = "blue", lwd = 2, lty = 2)
lines(grid, ci$up_ecdf, col = "blue", lwd = 2, lty = 2)
legend("topleft", col = c("black", "blue", "red"), lwd = 2, lty = c(1, 2, 2), legend = 
c("log-concave NPMLE", "CI for type = nrd", "CI for type = ECDFboot"), bty = "n")
}}
\keyword{htest}
\keyword{nonparametric}
