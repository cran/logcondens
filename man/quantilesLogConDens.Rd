\name{quantilesLogConDens}
\alias{quantilesLogConDens}
\title{Function to compute p0-Quantile of Fhat}
\description{Function to compute \eqn{p_0}-quantile of 

\deqn{\widehat F_m(t) = \int_{x_1}^t \widehat f_m(t) dt,}

where \eqn{\widehat f_m} is the log-concave density estimator, received via \code{\link{activeSetLogCon}}. The
formula to compute a quantile at \eqn{u \in [\widehat F_m(x_j), \widehat F_m(x_{j+1})]} for 
\eqn{j = 1, \ldots, n-1} is:

\deqn{\widehat F_m^{-1}(u) = x_j + (x_{j+1}-x_j)  G^{-1}_{(x_{j+1}-x_j)(\widehat \varphi_{j+1}-\widehat \varphi_j)} \Big( \frac{u - \widehat F_m(x_j)}{ \widehat F_m(x_{j+1}) -  \widehat F_m(x_j)}\Big),}

where \eqn{G^{-1}_\theta} is described in \code{\link{qloglin}}.}
\usage{quantilesLogConDens(p0, x, phi, Fhat)}
\arguments{
  \item{p0}{Real number where quantil should be computed.}
  \item{x}{Sorted vector of original observations \eqn{{\bold{x}} = (x_1, \ldots, x_m).}}
  \item{phi}{Vector \eqn{(\widehat \varphi_m(x_1), \ldots, \widehat \varphi_m(x_m)),} representing the function \eqn{\widehat \varphi_m}, as computed by \code{\link{activeSetLogCon}}.}
  \item{Fhat}{Vector \eqn{(\widehat F_{m,i})_{i=1}^m} with entries \deqn{\widehat F_{m,i} = \int_{x_1}^{x_i} \exp(\widehat \varphi_m(t)) dt,} as computed by \code{\link{activeSetLogCon}}.}
}
\value{Returns the real number \eqn{q_0 = \inf_{x}\{\widehat F_m(x) \ge p_0\}.}}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}

\examples{
## estimate gamma density
set.seed(1977)
x <- sort(rgamma(200, 2, 1))
res <- activeSetLogCon(x, w = NA, print = FALSE)

## compute 0.95 quantile of Fhat
q <- quantilesLogConDens(0.95, x, res$phi, res$Fhat)
plot(x, res$Fhat, type = 'l'); rug(x)
abline(h = 0.95, lty = 3); abline(v = q, lty = 3)
}
\keyword{htest}
\keyword{nonparametric}
