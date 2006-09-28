\name{quantilesLogConDens}
\alias{quantilesLogConDens}
\title{Function to compute p0-Quantile of F}
\description{Function to compute \eqn{p_0}-quantile of 

\deqn{\widehat F_m(t) = \int_{x_1}^t \widehat f_m(t) dt,}

where \eqn{\widehat f_m} is the log-concave density estimator, received via \code{\link{activeSetLogCon}}.}
\usage{quantilesLogConDens(p0, x, f, F, IsKnot)}
\arguments{
  \item{p0}{Real number where quantil should be computed.}
  \item{x}{Sorted vector of original observations \eqn{x = (x_1, \ldots, x_m).}}
  \item{f}{Vector \eqn{(\widehat f_m(x_1), \ldots, \widehat f_m(x_m)),} representing the function \eqn{\widehat f_m}, as computed by \code{\link{activeSetLogCon}}.}
  \item{F}{Vector \eqn{(\widehat F_{m,i})_{i=1}^m} with entries \deqn{\widehat F_{m,i} = \int_{x_1}^{x_i} \exp(\widehat \varphi_m(t)) dt,} as computed by \code{\link{activeSetLogCon}}.}
  \item{IsKnot}{Column vector with entries \cr
  IsKnot\eqn{_i = 1\{\widehat \varphi_m} has a kink at \eqn{x_i\}}, as computed by \code{\link{activeSetLogCon}}.}
}
\value{Returns the real number \eqn{q_0 = \inf_{x}\{\widehat F_m(x) \ge p_0\}.}}
\note{Since the log-density is piecewise linear, the corresponding distribution function can be analytically
computed from \eqn{x} and \eqn{\widehat \varphi_m}. The function \code{\link{quantilesLogConDens}} simply inverts this analytical function. However, for extreme quantiles
 (that are of course still smaller than \eqn{x_m}) and \eqn{m} big, this may be numerically unstable and should be 
 replaced by numerical search (e.g. bisection).}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}

\examples{
## estimate gamma density
set.seed(1977)
x <- sort(rgamma(200, 2, 1))
res <- activeSetLogCon(x, w = NA, print = FALSE)

## compute 0.95 quantile of estimated F
quantilesLogConDens(0.95, x, exp(res$phi), res$F, res$IsKnot)
}
\keyword{htest}
\keyword{nonparametric}
