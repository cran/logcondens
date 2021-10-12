\name{quantilesLogConDens}
\alias{quantilesLogConDens}
\title{Function to compute Quantiles of Fhat}
\description{Function to compute \eqn{p_0}-quantile of 

\deqn{\widehat F_m(t) = \int_{x_1}^t \widehat f_m(t) dt,}{\hat F_m(t) = \int_{x_1}^t \hat f_m(t) dt,}

where \eqn{\widehat f_m} is the log-concave density estimator, typically computed via \code{\link{logConDens}}
and \eqn{p_0} runs through the vector \code{ps}. 
The formula to compute a quantile at \eqn{u \in [\widehat F_m(x_j), \widehat F_m(x_{j+1})]}{u \in [\hat F_m(x_j), \hat F_m(x_{j+1})]} for 
\eqn{j = 1, \ldots, n-1} is:

\deqn{\widehat F_m^{-1}(u) = x_j + (x_{j+1}-x_j)  G^{-1}_{(x_{j+1}-x_j)(\widehat \phi_{j+1}-\widehat \phi_j)} \Big( \frac{u - \widehat F_m(x_j)}{ \widehat F_m(x_{j+1}) -  \widehat F_m(x_j)}\Big),}{\hat F_m^{-1}(u) = x_j + (x_{j+1}-x_j)  G^{-1}_{(x_{j+1}-x_j)(\hat \phi_{j+1}-\hat \phi_j)} ((u - \hat F_m(x_j))/(\hat F_m(x_{j+1}) - \hat F_m(x_j))),}

where \eqn{G^{-1}_\theta} is described in \code{\link{qloglin}}.}
\usage{quantilesLogConDens(ps, res)}
\arguments{
  \item{ps}{Vector of real numbers where quantiles should be computed.}
  \item{res}{An object of class \code{"dlc"}, usually a result of a call to \code{logConDens}.}
}
\value{Returns a data.frame with row \eqn{(p_{0, i}, q_{0, i})} where
\eqn{q_{0, i} = \inf_{x}\{\widehat F_m(x) \ge p_{0, i}\}}{q_{0, i} = inf_{x}\{\hat F_m(x) \ge p_{0, i}\}} and \eqn{p_{0, i}} runs through \code{ps}.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{https://www.imsv.unibe.ch/about_us/staff/prof_dr_duembgen_lutz/index_eng.html}}

\examples{
## estimate gamma density
set.seed(1977)
x <- rgamma(200, 2, 1)
res <- logConDens(x, smoothed = FALSE, print = FALSE)

## compute 0.95 quantile of Fhat
q <- quantilesLogConDens(0.95, res)[, "quantile"]
plot(res, which = "CDF", legend.pos = "none")
abline(h = 0.95, lty = 3); abline(v = q, lty = 3)
}
\keyword{htest}
\keyword{nonparametric}
